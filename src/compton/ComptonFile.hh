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

#include <string>
#include <vector>
#include <fstream>

namespace rtt_compton {

class ComptonFile {

  public:
  ComptonFile(std::string&, bool binary_=true);
  ~ComptonFile();

  // interface to data read:
  std::vector<std::vector<std::vector<std::vector<double>>>> 
  read_csk_data();

  // Accessors for evaluation points:
  std::vector<double> get_xi_pts() { return xi_pts; };
  std::vector<double> get_gout_pts() { return gout_pts; };
  std::vector<double> get_gin_pts() { return gin_pts; };
  std::vector<double> get_etemp_pts() { return etemp_pts; };

  private:
  // name of the file to be read
  std::string libfile;

  // input file stream:
  std::ifstream csk_data;

  // flag for binary file (default=true):
  bool binary;

  // evaluation points for the CSK (all except etemp read from data library):
  std::vector<double> gin_pts, gout_pts, xi_pts, etemp_pts;

  // read the csk library in binary format
  std::vector<std::vector<std::vector<std::vector<double>>>> 
  read_binary_csk_data();

  // read the csk library in ascii format
  std::vector<std::vector<std::vector<std::vector<double>>>> 
  read_ascii_csk_data();

};

} // end namespace rtt_compton

#endif
