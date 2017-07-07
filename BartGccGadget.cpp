/*******************************************************************
 * Description: Geometric Coil Compression (GCC) Gadget 
 * using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Author: Jia Sen
 * Lang: C++
 * Date: 05/17/2017
 * Version: 0.0.1
 *******************************************************************/

#include "BartGccGadget.h"

namespace Gadgetron {
  
  BartGccGadget::BartGccGadget() :
  image_counter_(0)
  {}

  int BartGccGadget::process_config(ACE_Message_Block* mb)
  {
    GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);
    
    // -------------------------------------------------
    
    ISMRMRD::IsmrmrdHeader h;
    try
    {
      deserialize(mb->rd_ptr(), h);
    }
    catch (...)
    {
      GDEBUG("Error parsing ISMRMRD Header");
    }
    
    size_t NE = h.encoding.size();
    num_encoding_spaces_ = NE;
    GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);
    
    
    return GADGET_OK;
  }
  
  int BartGccGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
  {
    
    if (perform_timing.value()) { gt_timer_local_.start("BartGccGadget::process"); }
    
    process_called_times_++;
    
    IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
    if (recon_bit_->rbit_.size() > num_encoding_spaces_)
    {
      GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
    }
    
    // for every encoding space
    for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
    {
      
      GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
      GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");
      
      
      // Check status of the folder containing the generated files (*.hdr & *.cfl)
      if (BartWorkingDirectory.value().empty()) {
	workLocation_ = workingDirectory.value();
      } else {
	workLocation_ = BartWorkingDirectory.value();
      }
      
      if (workLocation_.empty()) {
	GERROR("Undefined work location, bailing out\n");
	return GADGET_FAIL;
      }
      
      
      static std::string outputFolderPath = Gadgetron::CreateBartFileFolder(workLocation_);
      std::string generatedFilesFolder = std::string(outputFolderPath);
      boost::filesystem::path dir(generatedFilesFolder);
      if (boost::filesystem::create_directory(dir))
	GDEBUG("Folder to store *.hdr & *.cfl files is %s\n", generatedFilesFolder.c_str());
      else {
	GERROR("Folder to store *.hdr & *.cfl files doesn't exist...\n");
	return GADGET_FAIL;
      }
      
      //-------------------------------------------------------//
      // WRITE REFERENCE AND RAW DATA TO FILES 
      std::vector<uint16_t> DIMS_ref, DIMS;
      
      // Recon bit 
      Gadgetron::IsmrmrdReconBit it = recon_bit_->rbit_[e];
      
      // Grab a reference to the buffer containing the reference data
      auto  & dbuff_ref = it.ref_;
      hoNDArray< std::complex<float> >& ref = (*dbuff_ref).data_;
      
      // Grab a reference to the buffer containing the image data
      IsmrmrdDataBuffered & dbuff = it.data_;
      
      // kspace Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
      uint16_t E0 = static_cast<uint16_t>(dbuff.data_.get_size(0));
      uint16_t E1 = static_cast<uint16_t>(dbuff.data_.get_size(1));
      uint16_t E2 = static_cast<uint16_t>(dbuff.data_.get_size(2));
      uint16_t CHA = static_cast<uint16_t>(dbuff.data_.get_size(3));
      uint16_t N = static_cast<uint16_t>(dbuff.data_.get_size(4));
      uint16_t S = static_cast<uint16_t>(dbuff.data_.get_size(5));
      uint16_t LOC = static_cast<uint16_t>(dbuff.data_.get_size(6));
      DIMS = { E0, E1, E2, CHA, N, S, LOC };
      
      // prepare ref data for coil map calculation
      std::vector<size_t> data_dim;
      dbuff.data_.get_dimensions(data_dim);
      
      // reference Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
      uint16_t E0_ref = static_cast<uint16_t>(ref.get_size(0));
      uint16_t E1_ref = static_cast<uint16_t>(ref.get_size(1));
      uint16_t E2_ref = static_cast<uint16_t>(ref.get_size(2));
      uint16_t CHA_ref = static_cast<uint16_t>(ref.get_size(3));
      uint16_t N_ref = static_cast<uint16_t>(ref.get_size(4));
      uint16_t S_ref = static_cast<uint16_t>(ref.get_size(5));
      uint16_t LOC_ref = static_cast<uint16_t>(ref.get_size(6));
      DIMS_ref = { E0_ref, E1_ref, E2_ref, CHA_ref, N_ref, S_ref, LOC_ref };
      
      // Write reference data to disk
      GDEBUG_CONDITION_STREAM(true, "Reference Array [E0, E1, E2, CHA, N, S, LOC] = [" << E0_ref <<","<<E1_ref<<","<<E2_ref<<","<<CHA_ref<<","<<N_ref<<","<<S_ref<<","<<LOC_ref<<"]");
      write_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + "reference_data").c_str(),&ref);
      
      // Write kspace data to disk
      GDEBUG_CONDITION_STREAM(true, "Data Array [E0, E1, E2, CHA, N, S, LOC] = [" << E0 <<","<<E1<<","<<E2<<","<<CHA<<","<<N<<","<<S<<","<<LOC<<"]");
      write_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + "input_data").c_str(), &dbuff.data_);

      if ( ( DstChaNum.value() < CHA_ref) && (DstChaNum.value() < CHA) )
      {
	
	GDEBUG("Dst Channel Number is %d < %d \n", DstChaNum.value(), CHA);
	GDEBUG("Bart Geometric Coil Compression will be performed \n");
	//--------------------------------------------------------------------------//
	// Calling Bart Geometric Coil Compression
	std::ostringstream cmd2, cmd3, cmd4;
	std::replace(generatedFilesFolder.begin(), generatedFilesFolder.end(), '\\', '/');
	
	cmd2 << "/home/amax/bart/bart cc -r " << std::min<int>(CalibSize.value(), std::min<int>(E1_ref,E2_ref) )  << " -G " << " reference_data cc_matrix";	 
	GDEBUG("%s\n", cmd2.str().c_str());
	if (system(std::string("cd " + generatedFilesFolder + "&&" + cmd2.str()).c_str()))
	  return GADGET_FAIL;
	
	cmd3 << "/home/amax/bart/bart ccapply -p " << DstChaNum.value()  << " -G " << " input_data cc_matrix cc_input_data";	
	GDEBUG("%s\n", cmd3.str().c_str());
	if (system(std::string("cd " + generatedFilesFolder + "&&" + cmd3.str()).c_str()))
	  return GADGET_FAIL;
	
	cmd4 << "/home/amax/bart/bart ccapply -p " << DstChaNum.value()  << " -G " << " reference_data cc_matrix cc_reference_data";
	GDEBUG("%s\n", cmd4.str().c_str());
	if (system(std::string("cd " + generatedFilesFolder + "&&" + cmd4.str()).c_str()))
	  return GADGET_FAIL;
	
	//-------------------------------------------------------------------------//
	// Reformat the data back to gadgetron format by Bart command 
	std::string outputFile = "cc_input_data";
	std::string outputFile_ref = "cc_reference_data";
	
	boost::shared_ptr< hoNDArray<std::complex< float > > > DATA = 
	read_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + outputFile).c_str());
	
	boost::shared_ptr< hoNDArray<std::complex< float > > > REF = 
	read_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + outputFile_ref).c_str());
	
	// Delete Bart working BartWorkingDirectory
	if (BartWorkingDirectoryDelete.value()){
	  Gadgetron::cleanup(outputFolderPath);
	}
	
	//---------------------------------------------------------------------//
	// Assignment and pass down m1
	// Be careful about the memory 
	m1->getObjectPtr()->rbit_[e].data_.data_.clear();
	m1->getObjectPtr()->rbit_[e].ref_->data_.clear();
	
	m1->getObjectPtr()->rbit_[e].data_.data_ = *DATA.get();
	m1->getObjectPtr()->rbit_[e].ref_->data_ = *REF.get();
	
      }
      else{
	GDEBUG("Bart Geometric Coil Compression will be skipped \n");
      }
    

    }


    if (this->next()->putq(m1) < 0)
    {
      GERROR_STREAM("Put IsmrmrdReconData to Q failed ... ");
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(BartGccGadget)
}
