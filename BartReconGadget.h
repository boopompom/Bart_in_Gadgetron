/****************************************************************************************************************************
 * Description: Cartesian 3D Recon Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Recon chain: cc -g -> ecalib -> pics -> fakekspace -> fft -> image
 * Author: Jia Sen
 * Lang: C++
 * Date: 05/19/2017
 * Version: 0.0.1
 ****************************************************************************************************************************/

#ifndef BART_RECON_GADGET_H
#define BART_RECON_GADGET_H

#include "GenericReconGadget.h"
#include "GenericReconCartesianSpiritGadget.h"
#include "gadgetron_mricore_export.h"
#include "mri_core_data.h"
#include <gadgetron_paths.h>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include <utility>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include <iomanip>




#if defined (WIN32)
#ifdef __BUILD_GADGETRON_BartReconGadget__
#define EXPORTGADGETS_BartReconGadget __declspec(dllexport)
#else
#define EXPORTGADGETS_BartReconGadget __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_BartReconGadget
#endif


namespace Gadgetron {
  
  class EXPORTGADGETS_BartReconGadget BartReconGadget : public GenericReconGadget
  {
    
  public:
    GADGET_DECLARE(BartReconGadget);
    
    typedef GenericReconGadget BaseClass;
    typedef Gadgetron::GenericReconCartesianSpiritObj< std::complex<float> > ReconObjType;
    
    BartReconGadget();
    virtual ~BartReconGadget() = default;
    
  protected:
    GADGET_PROPERTY(BartWorkingDirectory, std::string, "Absolute path to temporary file location (will default to workingDirectory)", "/home/amax/");
    GADGET_PROPERTY(AbsoluteBartCommandScript_path, std::string, "Absolute path to bart script(s)", get_gadgetron_home() + "/share/gadgetron/bart");
    GADGET_PROPERTY(BartCommandScript_name, std::string, "Script file containing bart command(s) to be loaded", "");
    GADGET_PROPERTY(BartWorkingDirectoryDelete, bool, "Whether to delete BartWorkingDirectory", true);
    
    GADGET_PROPERTY(esp_map, int, "esp_map",2);
    GADGET_PROPERTY(n_iter_l1, int, "n_iter_l1", 15);
    GADGET_PROPERTY(lambda_l1, float, "lambda_l1", 0.002);
    
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
    
    long long image_counter_;
    std::string workLocation_;
    // record the recon kernel, coil maps etc. for every encoding space
    std::vector< ReconObjType > recon_obj_;
    
    // utility functions for read/write
    std::string CreateBartFileFolder(std::string workLocation_);
    template<typename U> void write_BART_Array(const char* filename, hoNDArray<U> *a);
    template <class T> boost::shared_ptr< hoNDArray<T> > read_BART_Array(const char* filename);
    std::string & getOutputFilename(const std::string & bartCommandLine);
    void cleanup(std::string &createdFiles);
    
    void perform_complex_coil_combine(ReconObjType& recon_obj);
    
    bool check_sampling_pattern(hoNDArray< std::complex< float > >  &out_data);

  };
  
  std::string BartReconGadget::CreateBartFileFolder(std::string workLocation_)
  {

    std::string outputFolderPath;
    
    time_t rawtime;
    char buff[80];
    time(&rawtime);
    strftime(buff, sizeof(buff), "%H_%M_%S__", localtime(&rawtime));
    std::mt19937::result_type seed = static_cast<unsigned long>(time(0));
    auto dice_rand = std::bind(std::uniform_int_distribution<int>(1, 10000), std::mt19937(seed));
    std::string time_id(buff + std::to_string(dice_rand()));
    outputFolderPath = workLocation_ + "bart_" + time_id + "/";

    return outputFolderPath;
  }
  
  template<typename U> void BartReconGadget::write_BART_Array(const char* filename, hoNDArray<U> *a)
  {
    std::vector<size_t> DIMS;
    for (int i = 0; i < a->get_number_of_dimensions(); i++)
    {
      DIMS.push_back(static_cast<size_t>(a->get_size(i)));
    }
    
    const size_t MAX_DIMS = 16;  // BART dims 16
    std::string filename_hdr = std::string(filename) + std::string(".hdr");
    std::vector<size_t> v(MAX_DIMS, 1);
    assert(DIMS.size() < MAX_DIMS);
    std::copy(DIMS.cbegin(), DIMS.cend(), v.begin());
    
    std::ofstream pFile_hdr;
    pFile_hdr.open(filename_hdr, std::ofstream::out);
    if (!pFile_hdr.is_open())
      GERROR("Failed to write into file: %s\n", filename);
    pFile_hdr << "# Dimensions\n";
    std::copy(v.cbegin(), v.cend(), std::ostream_iterator<size_t>(pFile_hdr, " "));
    pFile_hdr.close();
    
    
    std::string filename_s = std::string(filename) + std::string(".cfl");
    std::fstream pFile(filename_s, std::ios::out | std::ios::binary);
    if (!pFile.is_open())
      GERROR("Failed to write into file: %s\n", filename);
    
    pFile.write(reinterpret_cast<char*>(a->get_data_ptr()), a->get_number_of_elements()*sizeof(U));
    pFile.close();
  }
  
  
  template <class T> boost::shared_ptr< hoNDArray<T> > BartReconGadget::read_BART_Array(const char* filename)
  {
    
    std::string filename_hdr = std::string(filename) + std::string(".hdr");
    std::fstream infile_hdr(filename_hdr,std::ios::in | std::ios::binary);
    
    if (!infile_hdr.is_open())
      GERROR("Failed to open file: %s\n", filename_hdr.c_str());
    
    std::vector<size_t> DIMS;
    if (infile_hdr.is_open())
    {
      std::vector<std::string> tokens;
      std::string line;
      
      while (std::getline(infile_hdr, line, '\n'))
      {
	tokens.push_back(line);
      }
      
      // Parse the dimensions
      const std::string s = tokens[1];
      std::stringstream ss(s);
      std::string items;
      while (getline(ss, items, ' ')) {
	DIMS.push_back(std::stoi(items, nullptr, 10));
      }
      infile_hdr.close();
    }
    
    // convert from BART data of 16 dims to Gadgetron data of 7 dims
    // BART      dim order: [RO, E1, E2, CHA, MAP, TE, COEFF, COEFF2, ITER, CShift, Time1, Time2, Level, Slice, Avg]
    // Gadgetron dim order: [RO, E1, E2, CHA, N, S, LOC]
    std::vector<size_t> DIMS_GT;   
    DIMS_GT.push_back(DIMS[0]);    // RO
    DIMS_GT.push_back(DIMS[1]);    // E1
    DIMS_GT.push_back(DIMS[2]);    // E2
    DIMS_GT.push_back(DIMS[3]);    // CHA
    size_t dims_left = 1;
    for (int iter = 4; iter < DIMS.size(); iter ++)
      dims_left = dims_left*DIMS[iter];
    DIMS_GT.push_back(dims_left);
    DIMS_GT.push_back(1);
    DIMS_GT.push_back(1);
    
     // Load the cfl file
    std::string filename_s = std::string(filename) + std::string(".cfl");
    std::fstream infile(filename_s, std::ios::in | std::ios::binary);
    if (!infile.is_open()){
      GERROR("Failed to open file: %s\n", filename_s.c_str());
      return boost::shared_ptr< hoNDArray<T> >();
    }
    
    boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>(&DIMS_GT) );
    infile.read(reinterpret_cast<char*>(out->get_data_ptr()),sizeof(T)*out->get_number_of_elements());
    
    return out;
  }
  
  std::string & BartReconGadget::getOutputFilename(const std::string & bartCommandLine)
  {
    static std::vector<std::string> outputFile;
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char> > tokens(bartCommandLine, sep);
    for (auto itr = tokens.begin(); itr != tokens.end(); ++itr)
      outputFile.push_back(*itr);
    return (outputFile.back());
  }
  
  void BartReconGadget::cleanup(std::string &createdFiles)
  {
    boost::filesystem::remove_all(createdFiles);
  }

}
#endif //BART_RECON_GADGET_H
