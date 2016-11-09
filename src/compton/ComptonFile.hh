//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/ComptonFile.hh
 * \author Kendra Keady
 * \date   Fri Oct 14 2016
 * \brief  Header file for ComptonFile class
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */

#ifndef __compton_ComptonFile_hh__
#define __compton_ComptonFile_hh__

#include "ds++/SP.hh"
#include <string>
#include <vector>
#include <fstream>

namespace rtt_compton {

// forward declare the compton data class:
class ComptonData;

class ComptonFile {
  // useful typedef:
  typedef rtt_dsxx::SP<ComptonData> SP_CompData;

  public:
  //! Constructor with default binary library type
  ComptonFile(std::string&, bool binary_=true);
  //! Destructor
  ~ComptonFile();

  //! Public interface to compton library reader
  SP_CompData read_csk_data();

  //! Public interface to compton library reader
  SP_CompData read_lagrange_csk_data();


  private:
  // name of the file to be read
  std::string libfile;

  // input file stream for data processing
  std::ifstream csk_data;

  // flag for binary file (defaults to true in constructor)
  bool binary;

  // read the csk library in binary format
  SP_CompData read_binary_csk_data();

  // read the csk library in ascii format
  SP_CompData read_ascii_csk_data();

  // read the lagrange csk library in ascii format
  SP_CompData read_lagrange_ascii_csk_data();

  // read the lagrange csk library in binary format
  SP_CompData read_lagrange_binary_csk_data();

};

} // end namespace rtt_compton

#endif
