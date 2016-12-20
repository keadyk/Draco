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
#include <fstream>
#include <string>
#include <vector>

namespace rtt_compton {

// forward declare the compton data class:
class ComptonData;

class ComptonFile {
  // useful typedef:
  typedef rtt_dsxx::SP<ComptonData> SP_CompData;

public:
  //! Constructor with default binary library type
  ComptonFile(std::string &, bool binary_ = true);
  //! Destructor
  ~ComptonFile();

  //! Public interface to compton library reader-- returns a smart pointer
  //! to the data container
  SP_CompData read_mg_data();

private:
  // name of the file to be read
  std::string libfile;

  // input file stream for data processing
  std::ifstream csk_data;

  // flag for binary file (defaults to true in constructor)
  bool binary;

  SP_CompData read_mg_library_ascii(std::string &);
  SP_CompData read_mg_library_binary(std::string &);
};

} // end namespace rtt_compton

#endif
