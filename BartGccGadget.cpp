/****************************************************************************************************************************
 * Description: Geometric Coil Compression (GCC) Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Author: Jia Sen
 * Lang: C++
 * Date: 05/17/2017
 * Version: 0.0.1
 ****************************************************************************************************************************/

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
      std::vector<float> Temp_ref;
      Temp_ref.reserve(2 * E0_ref*E1_ref*E2_ref*CHA_ref*N_ref*S_ref*LOC_ref);
      GDEBUG_CONDITION_STREAM(true, "Reference Array [E0, E1, E2, CHA, N, S, LOC] = [" << E0_ref <<","<<E1_ref<<","<<E2_ref<<","<<CHA_ref<<","<<N_ref<<","<<S_ref<<","<<LOC_ref<<"]");
      
      for (uint16_t loc = 0; loc < LOC_ref; ++loc) {
	for (uint16_t s = 0; s < S_ref; ++s) {
	  for (uint16_t n = 0; n < N_ref; ++n) {
	    
	    //Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
	    //Each chunk will be [E0,E1,E2,CHA] big
	    std::vector<size_t> chunk_dims(4);
	    chunk_dims[0] = E0_ref;
	    chunk_dims[1] = E1_ref;
	    chunk_dims[2] = E2_ref;
	    chunk_dims[3] = CHA_ref;
	    // hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &(*dbuff_ref).data_(0, 0, 0, 0, n, s, loc));
	    hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &ref(0, 0, 0, 0, n, s, loc));
	    
	    std::vector<size_t> new_chunk_dims(1);
	    new_chunk_dims[0] = E0_ref*E1_ref*E2_ref*CHA_ref;
	    chunk.reshape(new_chunk_dims);
	    // Fill BART container
	    for (auto e : chunk) {
	      Temp_ref.push_back(e.real());
	      Temp_ref.push_back(e.imag());
	    }
	  }
	}	
	write_BART_Files(std::string(generatedFilesFolder + "reference_data").c_str(), DIMS_ref, Temp_ref);
      }
      
      // Write kspace data to disk
      std::vector<float> Temp;
      GDEBUG_CONDITION_STREAM(true, "Data Array [E0, E1, E2, CHA, N, S, LOC] = [" << E0 <<","<<E1<<","<<E2<<","<<CHA<<","<<N<<","<<S<<","<<LOC<<"]");
      Temp.reserve( 2*E0*E1*E2*CHA*N*S*LOC);
   
      for (uint16_t loc = 0; loc < LOC; loc++) {
	for (uint16_t s = 0; s < S; s++) {
	  for (uint16_t n = 0; n < N; n++) {
	    
	    //Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
	    //Each chunk will be [E0,E1,E2,CHA] big
	    std::vector<size_t> chunk_dims(4);
	    chunk_dims[0] = E0;
	    chunk_dims[1] = E1;
	    chunk_dims[2] = E2;
	    chunk_dims[3] = CHA;
	    hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &dbuff.data_(0, 0, 0, 0, n, s, loc));
	    
	    std::vector<size_t> new_chunk_dims(1);
	    new_chunk_dims[0] = E0*E1*E2*CHA;
	    chunk.reshape(new_chunk_dims);
	    // Fill BART container
	    for (auto e : chunk) {
	      Temp.push_back(e.real());
	      Temp.push_back(e.imag());
	    }
	  }
	}
      }
      
      write_BART_Files(std::string(generatedFilesFolder + "input_data").c_str(), DIMS, Temp);
      
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
      auto header = read_BART_hdr(std::string(generatedFilesFolder + outputFile).c_str());
     
      // READ FROM BART FILES
      std::pair< std::vector<size_t>, std::vector<std::complex<float> > > BART_DATA = Gadgetron::read_BART_files(std::string(generatedFilesFolder + outputFile).c_str());
      
      std::string outputFile_ref = "cc_reference_data";
      auto header_ref = read_BART_hdr(std::string(generatedFilesFolder + outputFile_ref).c_str());
     
      // READ FROM BART FILES
      std::pair< std::vector<size_t>, std::vector<std::complex<float> > > BART_REF = Gadgetron::read_BART_files(std::string(generatedFilesFolder + outputFile_ref).c_str());
      
      // Delete Bart working BartWorkingDirectory
      if (BartWorkingDirectoryDelete.value()){
	Gadgetron::cleanup(outputFolderPath);
      }
      
      
      // Grab data from BART files
      std::vector<size_t> BART_DATA_dims(1);
      BART_DATA_dims[0] = std::accumulate(BART_DATA.first.begin(), BART_DATA.first.end(), 1, std::multiplies<size_t>());
      hoNDArray<std::complex<float>> DATA = hoNDArray<std::complex<float>>(BART_DATA_dims, &BART_DATA.second[0]);
      
      // The image array data will be [E0,E1,E2,1,N,S,LOC]
      std::vector<size_t> data_dims(7);
      data_dims[0] = BART_DATA.first[0];
      data_dims[1] = BART_DATA.first[1];
      data_dims[2] = BART_DATA.first[2];
      data_dims[3] = BART_DATA.first[3];
      data_dims[4] = BART_DATA.first[4];
      data_dims[5] = BART_DATA.first[5];
      data_dims[6] = BART_DATA.first[6];
   
      std::vector<size_t> BART_REF_dims(1);
      BART_REF_dims[0] = std::accumulate(BART_REF.first.begin(), BART_REF.first.end(), 1, std::multiplies<size_t>());
      hoNDArray<std::complex<float>> REF = hoNDArray<std::complex<float>>(BART_REF_dims, &BART_REF.second[0]);
      std::vector<size_t> ref_dims(7);
      ref_dims[0] = BART_REF.first[0];
      ref_dims[1] = BART_REF.first[1];
      ref_dims[2] = BART_REF.first[2];
      ref_dims[3] = BART_REF.first[3];
      ref_dims[4] = BART_REF.first[4];
      ref_dims[5] = BART_REF.first[5];
      ref_dims[6] = BART_REF.first[6];
   
      DATA.reshape(data_dims);
      REF.reshape(ref_dims);
      
      
      std::vector<std::complex<float> > DATA_Final;
      DATA_Final.reserve(std::accumulate(data_dims.begin(), data_dims.end(), 1, std::multiplies<size_t>()));
      GDEBUG_CONDITION_STREAM(true, "[E0, E1, E2, CHA, N, S, LOC] = [" << data_dims[0] <<","<<data_dims[1] <<","<<data_dims[2] <<","<<data_dims[3] <<","<<data_dims[4] <<","<<data_dims[5] <<","<<data_dims[6] <<"]");

      std::vector<std::complex<float> > REF_Final;
      REF_Final.reserve(std::accumulate(ref_dims.begin(), ref_dims.end(), 1, std::multiplies<size_t>()));
      GDEBUG_CONDITION_STREAM(true, "[E0, E1, E2, CHA, N, S, LOC] = [" << ref_dims[0] <<","<<ref_dims[1] <<","<<ref_dims[2] <<","<<ref_dims[3] <<","<<ref_dims[4] <<","<<ref_dims[5] <<","<<ref_dims[6] <<"]");
      
      for (uint16_t loc = 0; loc < data_dims[6]; ++loc) {
	for (uint16_t s = 0; s < data_dims[5]; ++s) {
	  for (uint16_t n = 0; n < data_dims[4]; n += header[4]) {
	    
	    //Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
	    //Each chunk will be [E0,E1,E2,CHA] big
	    std::vector<size_t> chunk_dims(4), Temp_one_1d(1);
	    chunk_dims[0] = data_dims[0];
	    chunk_dims[1] = data_dims[1];
	    chunk_dims[2] = data_dims[2];
	    chunk_dims[3] = data_dims[3];
	    hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &DATA(0, 0, 0, 0, n, s, loc));
	    
	    Temp_one_1d[0] = chunk_dims[0] * chunk_dims[1] * chunk_dims[2] * chunk_dims[3];
	    chunk.reshape(Temp_one_1d);
	    DATA_Final.insert(DATA_Final.end(), chunk.begin(), chunk.end());
	    
	    //Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
	    //Each chunk will be [E0,E1,E2,CHA] big
	    std::vector<size_t> chunk_ref_dims(4), Ref_one_1d(1);
	    chunk_ref_dims[0] = ref_dims[0];
	    chunk_ref_dims[1] = ref_dims[1];
	    chunk_ref_dims[2] = ref_dims[2];
	    chunk_ref_dims[3] = ref_dims[3];
	    hoNDArray<std::complex<float> > chunk_ref = hoNDArray<std::complex<float> >(chunk_ref_dims, &REF(0, 0, 0, 0, n, s, loc));
	    
	    Ref_one_1d[0] = chunk_ref_dims[0] * chunk_ref_dims[1] * chunk_ref_dims[2] * chunk_ref_dims[3];
	    chunk_ref.reshape(Ref_one_1d);
	    REF_Final.insert(REF_Final.end(), chunk_ref.begin(), chunk_ref.end());
	    
	  }
	}
      }
      
      hoNDArray<std::complex< float >> data_cc_final(data_dims);
      std::copy(DATA_Final.begin(), DATA_Final.end(), data_cc_final.begin());
      
      hoNDArray<std::complex< float >> reference_cc_final(ref_dims);
      std::copy(REF_Final.begin(), REF_Final.end(), reference_cc_final.begin());
      
      //---------------------------------------------------------------------//
      // Assignment and pass down m1
      // Be careful about the memory 
      m1->getObjectPtr()->rbit_[e].data_.data_.clear();
      m1->getObjectPtr()->rbit_[e].ref_->data_.clear();
      m1->getObjectPtr()->rbit_[e].data_.data_ = data_cc_final;
      m1->getObjectPtr()->rbit_[e].ref_->data_ = reference_cc_final;

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
