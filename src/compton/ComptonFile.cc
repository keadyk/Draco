//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/ComptonFile.cc
 * \author Kendra Keady
 * \date   Fri Oct 14 2016
 * \brief  Implementation file for ComptonFile class.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 */
#include "ComptonFile.hh"
#include "ComptonData.hh"
#include "ds++/Assert.hh"
#include "ds++/Endian.hh"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>

namespace rtt_compton {

ComptonFile::ComptonFile(std::string &file_, bool binary_)
    : libfile(file_), binary(binary_) {}

ComptonFile::~ComptonFile() {}

// Interface to lagrange csk data read:
rtt_dsxx::SP<ComptonData> ComptonFile::read_mg_data() {
  // return value
  SP_CompData Cdata;
  if (binary) {
    Cdata = read_mg_library_binary(libfile);
  } else {
    Cdata = read_mg_library_ascii(libfile);
  }

  return Cdata;
}

rtt_dsxx::SP<ComptonData>
ComptonFile::read_mg_library_binary(std::string &infile) {
  {
    std::ostringstream msg;
    std::cout << "Binary read not implemented yet!" << std::endl;
    exit(1);
  }
}

rtt_dsxx::SP<ComptonData>
ComptonFile::read_mg_library_ascii(std::string &infile) {
  std::cout << "Reading multigroup CSK library from ascii file: " << infile
            << std::endl;

  // Open the output file:
  std::ifstream ascii_lib;
  ascii_lib.open(infile);

  // check file for "openness":
  if (!ascii_lib.is_open()) {
    std::cout << "Unable to open the specified input file!" << std::endl;
  }

  // FIRST read the number of electron temp break points and sub points,
  // the number of groups, and the number of legendre moments:
  std::string line;
  getline(ascii_lib, line);
  std::stringstream data_sizes(line);

  size_t netempbp;
  data_sizes >> netempbp;
  if (data_sizes.fail()) {
    Insist(0, "Failed read!");
  }
  std::cout << "number of etemp bps = " << netempbp << std::endl;

  size_t netemp;
  data_sizes >> netemp;
  if (data_sizes.fail()) {
    Insist(0, "Failed read!");
  }
  std::cout << "number of etemp evals = " << netemp << std::endl;

  size_t ngrp;
  data_sizes >> ngrp;
  if (data_sizes.fail()) {
    Insist(0, "Failed read!");
  }
  std::cout << "number of freq groups = " << ngrp << std::endl;

  size_t nleg;
  data_sizes >> nleg;
  if (data_sizes.fail()) {
    Insist(0, "Failed read!");
  }
  std::cout << "number of leg moments = " << nleg << std::endl;

  // initialize data container for the "raw" csk values:
  SP_CompData Cdata(new ComptonData(netemp, ngrp, nleg));

  // next line contains a list of electron temp breakpoints
  double etemp_bp;
  std::vector<double> etempbps;
  getline(ascii_lib, line);
  std::stringstream etemp_bps(line);
  for (size_t a = 0; a < netempbp; a++) {
    etemp_bps >> etemp_bp;
    if (etemp_bps.fail()) {
      Insist(0, "Failed read!");
    }
    etempbps.push_back(etemp_bp);
  }

  // next line contains a list of frequency group boundaries
  double grpbd;
  std::vector<double> grpbds;
  getline(ascii_lib, line);
  std::stringstream grp_bds(line);
  for (size_t a = 0; a <= ngrp; a++) {
    grp_bds >> grpbd;
    if (grp_bds.fail()) {
      Insist(0, "Failed read!");
    }
    grpbds.push_back(grpbd);
  }

  // stash these eval points in the comptondata container:
  Cdata->set_breakpts(etempbps);

  // make a temporary container for the raw data:
  std::vector<std::vector<std::vector<std::vector<double>>>> mg_csk_data(
      netemp, std::vector<std::vector<std::vector<double>>>(
                  ngrp, std::vector<std::vector<double>>(
                            ngrp, std::vector<double>(nleg, 0.0))));

  std::vector<double> etemppts;
  // now, for each electron temp eval point...
  for (size_t b = 0; b < netemp; b++) {
    // get the electron temp eval point:
    double etemp_pt;
    getline(ascii_lib, line);
    std::stringstream etemp_pts(line);
    etemp_pts >> etemp_pt;
    if (etemp_pts.fail()) {
      Insist(0, "Failed read!");
    }
    etemppts.push_back(etemp_pt);

    // loop over incident groups:
    for (size_t c = 0; c < ngrp; c++) {
      // loop over outgoing groups:
      for (size_t d = 0; d < ngrp; d++) {
        // get the group indices (and ignore them... we don't actually care)
        size_t gin, gout;
        getline(ascii_lib, line);
        std::stringstream data_line(line);
        data_line >> gin;
        data_line >> gout;
        std::cout << "gin and out " << gin << " " << gout << std::endl;
        // write all of the legendre moments for this g in/g out/etemp
        for (size_t e = 0; e < nleg; e++) {
          double data_point;
          data_line >> data_point;
          if (data_line.fail()) {
            Insist(0, "Failed read!");
          }
          mg_csk_data[b][c][d][e] = data_point;
          std::cout << "data point: " << data_point << std::endl;
        }
      }
    }
    // get a line and then toss it (blank)
    getline(ascii_lib, line);
  }

  Cdata->set_evalpts(etemppts, grpbds);
  Cdata->set_csk_data(mg_csk_data);
  // close the output file
  ascii_lib.close();

  // return the data smartpointer
  return Cdata;
}
}
