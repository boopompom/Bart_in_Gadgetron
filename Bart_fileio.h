#ifndef BART_FILEIO_H
#define BART_FILEIO_H

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#include <vector>
#include <cassert>
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

#pragma once

namespace Gadgetron{
  
  std::string CreateBartFileFolder(std::string workLocation_)
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
  
  template<typename T>
  void write_BART_hdr(const char* filename, std::vector<T> &DIMS)
  {
    const size_t MAX_DIMS = 16;
    std::string filename_s = std::string(filename) + std::string(".hdr");
    std::vector<size_t> v(MAX_DIMS, 1);
    assert(DIMS.size() < MAX_DIMS);
    std::copy(DIMS.cbegin(), DIMS.cend(), v.begin());
    std::ofstream pFile;
    pFile.open(filename_s, std::ofstream::out);
    if (!pFile.is_open())
      GERROR("Failed to write into file: %s\n", filename);
    pFile << "# Dimensions\n";
    std::copy(v.cbegin(), v.cend(), std::ostream_iterator<size_t>(pFile, " "));
    pFile.close();
  }
  
  template<typename T, typename U>
  void write_BART_Files(const char* filename, std::vector<T> &DIMS, std::vector<U> &DATA)
  {
    write_BART_hdr(filename, DIMS);
    std::string filename_s = std::string(filename) + std::string(".cfl");
    std::ofstream pFile(filename_s, std::ofstream::out | std::ofstream::binary);
    if (!pFile.is_open())
      GERROR("Failed to write into file: %s\n", filename);
    pFile.write(reinterpret_cast<char*>(&DATA[0]), DATA.size()*sizeof(float));
    pFile.close();
  }
  
  std::vector<size_t> read_BART_hdr(const char *filename)
  {
    std::string filename_s = std::string(filename) + std::string(".hdr");
    std::ifstream infile(filename_s);
    if (!infile.is_open())
      GERROR("Failed to open file: %s\n", filename_s.c_str());
    
    std::vector<size_t> DIMS;
    if (infile.is_open())
    {
      std::vector<std::string> tokens;
      std::string line;
      
      while (std::getline(infile, line, '\n'))
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
      infile.close();
    }
    
    return(DIMS);
  }
  
  
  std::pair< std::vector<size_t>, std::vector<std::complex<float> > >
  read_BART_files(const char *filename)
  {
    // Load the header file
    auto DIMS = read_BART_hdr(filename);
    
    // Load the cfl file
    std::string filename_s = std::string(filename) + std::string(".cfl");
    std::ifstream infile(filename_s, std::ifstream::binary);
    if (!infile.is_open())
      GERROR("Failed to open file: %s\n", filename_s.c_str());
    
    infile.seekg(0, infile.end);
    size_t size = static_cast<long>(infile.tellg());
    infile.seekg(0);
    std::vector<float> buffer(size);
    infile.read(reinterpret_cast<char*>(&buffer[0]), size);
    infile.close();
    // Reformat the data
    unsigned long count = 0;
    std::vector< std::complex<float> > Data;
    Data.reserve(buffer.size() / 2);
    
    for (size_t count = 0; count < buffer.size(); count += 2)
      Data.push_back(std::complex<float>(buffer[count], buffer[count + 1]));
    
    return (std::pair< std::vector<size_t>, std::vector< std::complex<float> > >(DIMS, Data));
  }
  
  std::string & getOutputFilename(const std::string & bartCommandLine)
  {
    static std::vector<std::string> outputFile;
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char> > tokens(bartCommandLine, sep);
    for (auto itr = tokens.begin(); itr != tokens.end(); ++itr)
      outputFile.push_back(*itr);
    return (outputFile.back());
  }
  
  void cleanup(std::string &createdFiles)
  {
    boost::filesystem::remove_all(createdFiles);
  }
  
}

#endif
