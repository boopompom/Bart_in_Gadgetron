/********************************************************************************************************************
 * Description: Cartesian 3D Recon Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Recon chain: cc -g -> ecalib -> pics -> fakekspace -> fft -> image
 * Author: Jia Sen
 * Lang: C++
 * Date: 05/19/2017
 * Version: 0.0.1
 ********************************************************************************************************************/

#include "BartReconGadget.h"
#include <armadillo>


namespace Gadgetron {
  
  BartReconGadget::BartReconGadget() : image_counter_(0)
  {}
  
  int BartReconGadget::process_config(ACE_Message_Block* mb)
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
    
    recon_obj_.resize(NE);
    
    
    return GADGET_OK;
  }
  
  int BartReconGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
  {
    
    if (perform_timing.value()) { gt_timer_local_.start("BartReconGadget::process"); }
    
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
      
      
      if (recon_bit_->rbit_[e].ref_)
      {
	
	// after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set	    
	if (perform_timing.value()) { gt_timer_.start("BartReconGadget::make_ref_coil_map"); }
	this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_,*recon_bit_->rbit_[e].data_.data_.get_dimensions(), recon_obj_[e].ref_calib_, recon_obj_[e].ref_coil_map_, e);
	if (perform_timing.value()) { gt_timer_.stop(); }
	
	// ----------------------------------------------------------
	
	// after this step, coil map is computed and stored in recon_obj_[e].coil_map_
	if (perform_timing.value()) { gt_timer_.start("BartReconGadget::perform_coil_map_estimation"); }
	this->perform_coil_map_estimation(recon_obj_[e].ref_coil_map_, recon_obj_[e].coil_map_, e);
	if (perform_timing.value()) { gt_timer_.stop(); }
	
      }
      
      //-------------------------Bart Recon Start-------------------------------------//
      // Check status of bart commands script 
      std::string CommandScript = AbsoluteBartCommandScript_path.value() + "/" + BartCommandScript_name.value();
      if (!boost::filesystem::exists(CommandScript))
      {
	GERROR("Can't find bart commands script: %s!\n", CommandScript.c_str());
	return GADGET_FAIL;
      }
      GDEBUG("Bart Command Script: %s\n", CommandScript.c_str());
      
      
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
      
      static std::string outputFolderPath = CreateBartFileFolder(workLocation_);
      std::string generatedFilesFolder = std::string(outputFolderPath);
      boost::filesystem::path dir(generatedFilesFolder);
      if (boost::filesystem::create_directory(dir))
	GDEBUG("Folder to store *.hdr & *.cfl files is %s\n", generatedFilesFolder.c_str());
      else {
	GERROR("Folder to store *.hdr & *.cfl files doesn't exist...\n");
	return GADGET_FAIL;
      }
      
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
      GDEBUG_CONDITION_STREAM(verbose.value(), "Reference Array [E0, E1, E2, CHA, N, S, LOC] = [" << E0_ref <<","<<E1_ref<<","<<E2_ref<<","<<CHA_ref<<","<<N_ref<<","<<S_ref<<","<<LOC_ref<<"]");
      write_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + "reference_data").c_str(),&ref);
      
      // Write kspace data to disk
      if (perform_timing.value()) { gt_timer_.start("BartReconGadget::write out kspace array to disk"); }
      
      GDEBUG_CONDITION_STREAM(verbose.value(), "Data Array [E0, E1, E2, CHA, N, S, LOC] = [" << E0 <<","<<E1<<","<<E2<<","<<CHA<<","<<N<<","<<S<<","<<LOC<<"]");
      write_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + "input_data").c_str(), &dbuff.data_);
      
      if (perform_timing.value()) { gt_timer_.stop(); } 
      
      
      // Run Bart script
      int isvdPDS = check_sampling_pattern(recon_bit_->rbit_[e].data_.data_);
      if (perform_timing.value()) { gt_timer_.start("BartReconGadget::ESPIRiT calibration + PICS reconstruction"); }
      std::string Commands_Line, last_command;
      ifstream inputFile;
      
      inputFile.open(CommandScript);
      if (inputFile.is_open())
      {
	while (getline(inputFile, Commands_Line))
	{
	  if(Commands_Line.empty() || Commands_Line.find("bart") == std::string::npos )
	    continue;
	  last_command = Commands_Line;
	}
	
	std::ostringstream Script_params;
	Script_params<<" -w "<< lambda_l1.value()<<" -i " << n_iter_l1.value() <<" -m "<< esp_map.value() <<" input_data "; 

	auto ret = system(std::string("cd " + generatedFilesFolder + "&&" + CommandScript + Script_params.str()).c_str()); 
        (void)ret;

	inputFile.close();
      }
      else
      {
	GERROR("Unable to open %s\n", CommandScript.c_str());
	cleanup(outputFolderPath);
	return GADGET_FAIL;
      }
      
      // Grab data from BART files
      std::string outputFile = getOutputFilename(last_command);
      
      boost::shared_ptr< hoNDArray<std::complex< float > > > DATA = 
      read_BART_Array<std::complex< float > >(std::string(generatedFilesFolder + outputFile).c_str());
      
      if (BartWorkingDirectoryDelete.value())
      {
	cleanup(outputFolderPath);
      }
      if (perform_timing.value()) { gt_timer_.stop(); } 
      //-------------------------Bart Recon Finished-------------------------------------//
      
      // Coil Combination
      recon_obj_[e].full_kspace_ = *DATA.get();
      if (this->perform_timing.value()) gt_timer_.start("BartReconGadget::perform_coil_combination using CSM... ");
      this->perform_complex_coil_combine(recon_obj_[e]);
      if (this->perform_timing.value()) gt_timer_.stop();
      
      // sending out image array
      if (recon_obj_[e].recon_res_.data_.get_number_of_elements() > 0)
      {
	
	if (perform_timing.value()) { gt_timer_.start("BartReconGadget::compute_image_header"); }
	this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e);
	if (perform_timing.value()) { gt_timer_.stop(); }
	
	
	if (perform_timing.value()) { gt_timer_.start("BartReconGadget::send_out_image_array"); }
	std::ostringstream ostr_image;
	ostr_image << "ESPIRIT-" << std::setprecision(5) << lambda_l1.value();
	std::string imageInfo = ostr_image.str();
	this->send_out_image_array_online(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR,imageInfo);
	if (perform_timing.value()) { gt_timer_.stop(); }
	
      }
      
      recon_bit_->rbit_[e].ref_ = boost::none;
      recon_obj_[e].recon_res_.data_.clear();
      recon_obj_[e].recon_res_.headers_.clear();
      recon_obj_[e].recon_res_.meta_.clear();
    }
    
    m1->release();
    return GADGET_OK;
  }
  
   bool BartReconGadget::check_sampling_pattern(hoNDArray< std::complex< float > >  &out_data)
  {
    size_t RO = out_data.get_size(0);
    size_t E1 = out_data.get_size(1);
    size_t E2 = out_data.get_size(2);
    size_t CHA = out_data.get_size(3);
    size_t N = out_data.get_size(4);
    size_t S = out_data.get_size(5);
    size_t SLC = out_data.get_size(6);
    
    double maxr = 0;
    arma::Mat<int> mask(E1,E2,arma::fill::zeros);
    
    size_t e1, e2, n, s;
    size_t num_readout_lines = 0;
    for (s = 0; s < S; s++)
    {
      for (n = 0; n < N; n++)
      {
	for (e2 = 0; e2 < E2; e2++)
	{
	  for (e1 = 0; e1 < E1; e1++)
	  {
	    if (std::abs(out_data(RO / 2, e1, e2, 0, n)) > 0)
	    {
	      num_readout_lines++;
	      mask(e1,e2) = 1;
	    }
	  }
	}
      }
    }
    
    bool is_ipat2D = 1;
    bool is_vdPDS = 0;
    if (num_readout_lines > 0)
    {
      float effective_acce_factor = float(S*N*E1*E2)/float(num_readout_lines);
      GDEBUG_STREAM("effective_acce_factor : " << effective_acce_factor);
      
      
      arma::uvec E1_center_acq = arma::find(mask.row(E1/2));
      arma::uvec E2_center_acq = arma::find(mask.col(E2/2));
      arma::vec E1_delta(E1_center_acq.n_elem - 1,arma::fill::zeros);
      arma::vec E2_delta(E2_center_acq.n_elem - 1,arma::fill::zeros);
      for (int iter = 0; iter < E1_center_acq.n_elem - 1; iter++)
      {
	E1_delta(iter) = E1_center_acq(iter + 1) - E1_center_acq(iter);
      }
      for (int iter = 0; iter < E2_center_acq.n_elem - 1; iter++)
      {
	E2_delta(iter) = E2_center_acq(iter + 1) - E2_center_acq(iter);
      }
      arma::uvec E1_delta_R = arma::find(E1_delta > 1.5*acceFactorE1_[0]);
      arma::uvec E2_delta_R = arma::find(E2_delta > 1.5*acceFactorE2_[0]);
      if ( E1_delta_R.n_elem > 5 || E2_delta_R.n_elem > 5 )
      {
	is_ipat2D = 0;
	is_vdPDS = 1;
      }
      GDEBUG_STREAM("Is variable density sampling : " << is_vdPDS);
      
    }
   
    return is_vdPDS;
    
  }
  
  void BartReconGadget::perform_complex_coil_combine(ReconObjType& recon_obj)
  {
    try
    {
      size_t RO = recon_obj.full_kspace_.get_size(0);
      size_t E1 = recon_obj.full_kspace_.get_size(1);
      size_t E2 = recon_obj.full_kspace_.get_size(2);
      size_t dstCHA = recon_obj.full_kspace_.get_size(3);
      size_t N = recon_obj.full_kspace_.get_size(4);
      size_t S = recon_obj.full_kspace_.get_size(5);
      size_t SLC = recon_obj.full_kspace_.get_size(6);
      
      recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);
      Gadgetron::clear(recon_obj.recon_res_.data_);
      
      if (E2>1)
      {
	Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_obj.full_kspace_, complex_im_recon_buf_);
      }
      else
      {
	Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_obj.full_kspace_, complex_im_recon_buf_);
      }
      
      size_t num = N*S*SLC;
      long long ii;
      
      #pragma omp parallel default(none) private(ii) shared(num, N, S, recon_obj, RO, E1, E2, dstCHA) if(num>1)
      {
	hoNDArray< std::complex<float> > complexImBuf(RO, E1, E2, dstCHA);
	
	#pragma omp for
	for (ii = 0; ii < num; ii++)
	{
	  size_t slc = ii / (N*S);
	  size_t s = (ii - slc*N*S) / N;
	  size_t n = ii - slc*N*S - s*N;
	  
	  size_t coilMapN = n;
	  if (coilMapN >= recon_obj.coil_map_.get_size(5)) coilMapN = recon_obj.coil_map_.get_size(5) - 1;
	  
	  size_t coilMapS = s;
	  if (coilMapS >= recon_obj.coil_map_.get_size(6)) coilMapS = recon_obj.coil_map_.get_size(6) - 1;
	  
	  hoNDArray< std::complex<float> > complexIm(RO, E1, E2, dstCHA, &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc)));
	  hoNDArray< std::complex<float> > coilMap(RO, E1, E2, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, coilMapN, coilMapS, slc)));
	  hoNDArray< std::complex<float> > combined(RO, E1, E2, 1, &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc)));
	  
	  Gadgetron::multiplyConj(complexIm, coilMap, complexImBuf);
	  Gadgetron::sum_over_dimension(complexImBuf, combined, 3);
	}
      }
    }
    catch (...)
    {
      GADGET_THROW("Errors happened in BartReconGadget::perform_complex_coil_combine(...) ... ");
    }
  }
  
  GADGET_FACTORY_DECLARE(BartReconGadget)
}
