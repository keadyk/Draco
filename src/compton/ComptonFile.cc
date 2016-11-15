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

// Interface to csk data read:
rtt_dsxx::SP<ComptonData> ComptonFile::read_csk_data() {
  // return value
  SP_CompData Cdata;
  if (binary) {
    Cdata = read_binary_csk_data();
  } else {
    Cdata = read_ascii_csk_data();
  }

  return Cdata;
}

// Interface to lagrange csk data read:
//std::vector<std::vector<std::vector<std::vector<double>>>>
rtt_dsxx::SP<ComptonData> ComptonFile::read_lagrange_csk_data() {
  // return value
  SP_CompData Cdata;
  if (binary) {
    Cdata = read_lagrange_binary_csk_data();
  } else {
    Cdata = read_lagrange_ascii_csk_data();
  }

  return Cdata;
}

// TODO: IF we go with the angular moment-based method, we should take out the
// xi-eval-pt read (line 2 of the file)
rtt_dsxx::SP<ComptonData> ComptonFile::read_binary_csk_data() {
  // *************************************************************************
  // ************************ OPEN THE DATA FILE *****************************
  // *************************************************************************
  // open the file:
  csk_data.open(libfile, std::ios::binary);

  // *************************************************************************
  // ********************* READ SIZES/ALLOCATE ARRAYS ************************
  // *************************************************************************
  // first line of the file tells us how to allocate the raw data array
  // (which is the return value of the function:)
  if (!csk_data.is_open()) {
    std::ostringstream msg;
    msg << "ACK, we failed to open " << libfile << "! ";
    Insist(0, msg.str());
  }

  // *************************************************************************
  // ************************* GET DATA ARRAY SIZES **************************
  // *************************************************************************
  // declare variables for number of eval pts in gamma in, gamma out, xi, and
  // etemp:
  size_t n_gin, n_gout, n_xi, n_etemp;

  // make a vector to hold the raw size data:
  std::vector<char> size_data(4 * sizeof(size_t));

  // grab the first value (number of electron temp points)
  csk_data.read(&size_data[0], 4 * sizeof(size_t));

  // check for failure
  if (csk_data.fail()) {
    Insist(0, "Failed to read CSK data sizes!");
  }

  // cast raw character data to size_t
  std::memcpy(&n_etemp, &size_data[0], sizeof(size_t));

  // cast raw character data to size_t
  std::memcpy(&n_gin, &size_data[sizeof(size_t)], sizeof(size_t));

  // cast raw character data to size_t
  std::memcpy(&n_gout, &size_data[2 * sizeof(size_t)], sizeof(size_t));

  // cast raw character data to size_t
  std::memcpy(&n_xi, &size_data[3 * sizeof(size_t)], sizeof(size_t));

  // initialize data container for the "raw" csk values:
  bool lagrange = false;
  SP_CompData Cdata(new ComptonData(n_etemp, n_gin, n_gout, n_xi, lagrange));

  // *************************************************************************
  // ************************ GET THE XI POINTS ******************************
  // *************************************************************************
  // get the next n_xi * sizeof(double) bytes from the file
  // allocate temporary vectors for xi data:
  std::vector<double> xi_pts(n_xi, 0.0);
  std::vector<char> xi_data(n_xi * sizeof(double));
  csk_data.read(&xi_data[0], n_xi * sizeof(double));
  // check for failure:
  if (csk_data.fail()) {
    Insist(0, "Failed to read xi eval points!");
  }

  for (size_t m = 0; m < n_xi; m++) {
    double xi;
    std::memcpy(&xi, &xi_data[m * sizeof(double)], sizeof(double));
    xi_pts[m] = xi;
  }

  // *************************************************************************
  // ************************** READ THE DATA ********************************
  // *************************************************************************
  //make temporary vectors for the evaluation points:
  std::vector<double> gin_pts(n_gin, 0.0);
  std::vector<double> gout_pts(n_gout, 0.0);
  std::vector<double> etemp_pts(n_etemp, 0.0);

  std::vector<std::vector<std::vector<std::vector<double>>>> raw_csk_data(
      n_etemp, std::vector<std::vector<std::vector<double>>>(
                   n_gin, std::vector<std::vector<double>>(
                              n_gout, std::vector<double>(n_xi, 0.0))));

  for (size_t a = 0; a < n_etemp; a++) // for each electron temp expected
  {
    // get the electron temp for this block:
    double etemp;
    std::vector<char> etemp_val(sizeof(double));
    csk_data.read(&etemp_val[0], sizeof(double));
    memcpy(&etemp, &etemp_val[0], sizeof(double));
    etemp_pts[a] = etemp;

    // get the REST of the data for the block:
    std::vector<char> data_block(n_gin * n_gout * (n_xi + 2) * sizeof(double));
    csk_data.read(&data_block[0], n_gin * n_gout * (n_xi + 2) * sizeof(double));

    for (size_t b = 0; b < n_gin; b++) // for each gamma in
    {
      for (size_t c = 0; c < n_gout; c++) // for each gamma out
      {
        // calculate the offset into the data array
        size_t offset = (b * n_gout + c) * (n_xi + 2) * sizeof(double);
        // get the first two values from the data block:
        double gin, gout;
        memcpy(&gin, &data_block[offset], sizeof(double));
        memcpy(&gout, &data_block[offset + sizeof(double)], sizeof(double));
        // store the gamma in points for the first electron temperature:
        if (a == 0) {
          if (c == 0) {
            gin_pts[b] = gin;
          }
          // store the gamma out points for the first gamma in:
          if (b == 0) {
            gout_pts[c] = gout;
          }
        }
        // update the data offset
        offset += 2 * sizeof(double);
        for (size_t d = 0; d < n_xi; d++) // for each xi value
        {
          double csk_value;
          // copy from the data array
          std::memcpy(&csk_value, &data_block[offset + d * sizeof(double)],
                      sizeof(double));
          raw_csk_data[a][b][c][d] = csk_value;

        } // end xi loop
      }   // end gamma out loop
    }     // end gamma in loop
  }       // end electron temperature loop

  // set the evaluation points now that we've read everything
  Cdata->set_evalpts(etemp_pts, gin_pts, gout_pts, xi_pts);
  Cdata->set_csk_data(raw_csk_data);
  // *************************************************************************
  // ************************** CLOSE THE FILE *******************************
  // *************************************************************************
  csk_data.close();

  std::cout << "CSK data read successfully!" << std::endl;
  return Cdata;
}

// read formatted csk data from a single ascii CSK library file
// this method assumes the data file is formatted in the following way:
// ----------------------------(Beginning of file)-----------------------------
// <# electron temp pts> <# gamma in pts> <# gamma out pts> <# xi pts> \n
// <xi 1> <xi 2> <xi 3 > ... <xi last>
// <electron temp 1>
// <gamma in 1> <gamma out 1> <CSK>
// <gamma in 1> <gamma out 2> <CSK>
// <gamma in 1> <gamma out 3> <CSK>
// ...
// <gamma in 1> <gamma out last> <CSK>
// <gamma in 2> <gamma out 1> <CSK>
// <gamma in 2> <gamma out 2> <CSK>
// <gamma in 2> <gamma out 3> <CSK>
// ...
// <gamma in 2> <gamma out last> <CSK>
// ...
// <gamma in last> <gamma out 1> <CSK>
// <gamma in last> <gamma out 2> <CSK>
// <gamma in last> <gamma out 3> <CSK>
// ...
// <gamma in last> <gamma out last> <CSK>
// <electron temp 2>
//
// -------------------------------(End of file)---------------------------------
// TODO: IF we go with the angular moment-based method, we should take out the
// xi-eval-pt read (line 2 of the file)
rtt_dsxx::SP<ComptonData> ComptonFile::read_ascii_csk_data() {
  // *************************************************************************
  // ************************* OPEN THE DATA FILE ****************************
  // *************************************************************************
  // open the file:
  csk_data.open(libfile);

  // *************************************************************************
  // ********************* READ SIZES/ALLOCATE ARRAYS ************************
  // *************************************************************************
  // first line of the file tells us how to allocate the raw data array
  // (which is the return value of the function:)
  if (!csk_data.is_open()) {
    std::ostringstream msg;
    msg << "ACK, we failed to open " << libfile << "! ";
    Insist(0, msg.str());
  }

  // *************************************************************************
  // ************************* GET DATA ARRAY SIZES **************************
  // *************************************************************************
  // declare variables for number of eval pts in gamma in, gamma out, xi, and
  // etemp:
  size_t n_gin, n_gout, n_xi, n_etemp;

  // string to hold the first line:
  std::string sizeline;

  // grab the first line from the datafile:
  getline(csk_data, sizeline);

  // stringstream to pipe in values from the line:
  std::stringstream size_data(sizeline);

  // grab the first value (number of electron temp points)
  size_data >> n_etemp;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of electron temperature points!");
  }

  // grab the second value (number of gamma ins)
  size_data >> n_gin;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of gamma in points!");
  }

  // grab the third value (number of gamma outs)
  size_data >> n_gout;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of gamma out points!");
  }

  // grab the fourth value (number of xi points)
  size_data >> n_xi;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of xi points!");
  }

  // initialize data container for the "raw" csk values:
  bool lagrange = false;
  SP_CompData Cdata(new ComptonData(n_etemp, n_gin, n_gout, n_xi, lagrange));

  // *************************************************************************
  // ************************ GET THE XI POINTS ******************************
  // *************************************************************************
  // make a temporary vector to hold the xi points
  std::vector<double> xi_pts(n_xi, 0.0);
  // (recycle the string from the previous line)
  getline(csk_data, sizeline);
  std::stringstream xi_data(sizeline);
  for (size_t m = 0; m < n_xi; m++) {
    double xi;
    xi_data >> xi;
    if (!xi_data.fail()) {
      xi_pts[m] = xi;
    } else {
      Insist(0, "CSK data read failed!");
    }
  }

  // *************************************************************************
  // ************************** READ THE DATA ********************************
  // *************************************************************************
  // declare a string representing a single line from the file:
  std::string file_line;

  //make temporary vectors for the evaluation points:
  std::vector<double> gin_pts(n_gin, 0.0);
  std::vector<double> gout_pts(n_gout, 0.0);
  std::vector<double> etemp_pts(n_etemp, 0.0);

  std::vector<std::vector<std::vector<std::vector<double>>>> raw_csk_data(
      n_etemp, std::vector<std::vector<std::vector<double>>>(
                   n_gin, std::vector<std::vector<double>>(
                              n_gout, std::vector<double>(n_xi, 0.0))));

  for (size_t a = 0; a < n_etemp; a++) // for each electron temp expected
  {
    double etemp;
    getline(csk_data, file_line);
    // declare a stringstream to extract data from the line:
    std::stringstream line_data(file_line);
    line_data >> etemp;
    etemp_pts[a] = etemp;
    for (size_t b = 0; b < n_gin; b++) // for each gamma in
    {
      for (size_t c = 0; c < n_gout; c++) // for each gamma out
      {
        // get a line of data from the file
        getline(csk_data, file_line);
        std::stringstream line_data(file_line);

        // get the first two values from the line:
        double gin, gout;
        line_data >> gin;
        line_data >> gout;

        // store the gamma in points for the first electron temperature:
        if (a == 0) {
          if (c == 0) {
            if (!line_data.fail()) {
              gin_pts[b] = gin;
            } else {
              Insist(0, "CSK data read failed!");
            }
          }
          // store the gamma out points for the first gamma in:
          if (b == 0) {
            if (!line_data.fail()) {
              gout_pts[c] = gout;
            } else {
              Insist(0, "CSK data read failed!");
            }
          }
        }

        for (size_t d = 0; d < n_xi; d++) // for each xi value
        {
          double csk_value;
          // get the data value and store it!
          line_data >> csk_value;
          if (!line_data.fail()) {
            raw_csk_data[a][b][c][d] = csk_value;
          } else {
            Insist(0, "CSK data read failed!");
          }

        } // end xi loop
      }   // end gamma out loop
    }     // end gamma in loop
  }       // end electron temperature loop

  Cdata->set_evalpts(etemp_pts, gin_pts, gout_pts, xi_pts);
  Cdata->set_csk_data(raw_csk_data);
  // *************************************************************************
  // ************************** CLOSE THE FILE *******************************
  // *************************************************************************
  csk_data.close();

  std::cout << "CSK data read successfully!" << std::endl;

  return Cdata;
}

// TODO: IF we go with the angular moment-based method, we should take out the
// xi-eval-pt read (line 2 of the file). Also, this routine assumes the etemp,
// gin, and gout variables will all be treated using lagrange interpolation
rtt_dsxx::SP<ComptonData> ComptonFile::read_lagrange_ascii_csk_data() {
  // *************************************************************************
  // ************************* OPEN THE DATA FILE ****************************
  // *************************************************************************
  // open the file:
  csk_data.open(libfile);

  // *************************************************************************
  // ********************* READ SIZES/ALLOCATE ARRAYS ************************
  // *************************************************************************
  // first line of the file tells us how to allocate the raw data array
  // (which is the return value of the function:)
  if (!csk_data.is_open()) {
    std::ostringstream msg;
    msg << "ACK, we failed to open " << libfile << "! ";
    Insist(0, msg.str());
  }

  // *************************************************************************
  // ************************* GET DATA ARRAY SIZES **************************
  // *************************************************************************
  // declare variables for number of eval pts in gamma in, gamma out, xi, and
  // etemp (including breakpoints and "local" eval points):
  size_t n_ginbp, n_goutbp, n_etempbp;
  size_t n_gin, n_gout, n_xi, n_etemp;

  // string to hold the first line:
  std::string sizeline;

  // grab the first line from the datafile:
  getline(csk_data, sizeline);

  // stringstream to pipe in values from the line:
  std::stringstream size_data(sizeline);

  // grab the first two values (number of electron break pts/local pts)
  size_data >> n_etempbp;
  size_data >> n_etemp;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of electron temperature points!");
  }

  // grab the second value (number of gamma ins)
  size_data >> n_ginbp;
  size_data >> n_gin;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of gamma in points!");
  }

  // grab the third value (number of gamma outs)
  size_data >> n_goutbp;
  size_data >> n_gout;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of gamma out points!");
  }

  // grab the fourth value (number of xi points)
  size_data >> n_xi;
  // check for failure
  if (size_data.fail()) {
    Insist(0, "Failed to read number of xi points!");
  }

  // initialize data container for the "raw" csk values:
  SP_CompData Cdata;
  Cdata.reset(new ComptonData(n_etemp, n_gin, n_gout, n_xi));

  // *************************************************************************
  // *********************** GET THE BREAKPOINTS *****************************
  // *************************************************************************

  // (recycle the string from the previous line)
  std::vector<double> etemp_breakpts(n_etempbp, 0.0);
  getline(csk_data, sizeline);
  std::stringstream etemp_data(sizeline);
  for (size_t m = 0; m < n_etempbp; m++) {
    double bp;
    etemp_data >> bp;
    if (!etemp_data.fail()) {
      etemp_breakpts[m] = bp;
    } else {
      Insist(0, "CSK data read failed!");
    }
  }

  std::vector<double> gin_breakpts(n_ginbp, 0.0);
  getline(csk_data, sizeline);
  std::stringstream gin_data(sizeline);
  for (size_t m = 0; m < n_ginbp; m++) {
    double bp;
    gin_data >> bp;
    if (!gin_data.fail()) {
      gin_breakpts[m] = bp;
    } else {
      Insist(0, "CSK data read failed!");
    }
  }

  std::vector<double> gout_breakpts(n_goutbp, 0.0);
  getline(csk_data, sizeline);
  std::stringstream gout_data(sizeline);
  for (size_t m = 0; m < n_goutbp; m++) {
    double bp;
    gout_data >> bp;
    if (!gout_data.fail()) {
      gout_breakpts[m] = bp;
    } else {
      Insist(0, "CSK data read failed!");
    }
  }

  // All breakpoints are now read-- set them in the Cdata object
  Cdata->set_breakpts(etemp_breakpts, gin_breakpts, gout_breakpts);

  // *************************************************************************
  // ************************ GET THE XI POINTS ******************************
  // *************************************************************************
  std::vector<double> xi_pts(n_xi, 0.0);
  // (recycle the string from the previous line)
  getline(csk_data, sizeline);
  std::stringstream xi_data(sizeline);
  for (size_t m = 0; m < n_xi; m++) {
    double xi;
    xi_data >> xi;
    if (!xi_data.fail()) {
      xi_pts[m] = xi;
    } else {
      Insist(0, "CSK data read failed!");
    }
  }

  // *************************************************************************
  // ************************** READ THE DATA ********************************
  // *************************************************************************
  // declare a string representing a single line from the file:
  std::string file_line;

  {
    std::vector<double> etemp_pts(n_etemp, 0.0);
    std::vector<double> gin_pts(n_gin, 0.0);
    std::vector<double> gout_pts(n_gout, 0.0);

    std::vector<std::vector<std::vector<std::vector<double>>>> raw_csk_data(
        n_etemp, std::vector<std::vector<std::vector<double>>>(
                     3 * n_gin, std::vector<std::vector<double>>(
                                    n_gout, std::vector<double>(n_xi, 0.0))));

    for (size_t a = 0; a < n_etemp; a++) // for each electron temp expected
    {
      double etemp;
      getline(csk_data, file_line);
      // declare a stringstream to extract data from the line:
      std::stringstream line_data(file_line);
      line_data >> etemp;
      etemp_pts[a] = etemp;
      // For each electron temperature, there are THREE regions of data,
      // based on the location of the boundary layers in gamma in/gamma out
      for (size_t k = 0; k < 3; k++) {
        for (size_t b = 0; b < n_gin; b++) // for each gamma in
        {
          for (size_t c = 0; c < n_gout; c++) // for each gamma out
          {
            // get a line of data from the file
            getline(csk_data, file_line);
            std::stringstream line_data(file_line);

            // get the first two values from the line:
            double gin, gout;
            line_data >> gin;
            line_data >> gout;
            // store the gamma in points for the first electron temperature:
            if (a == 0) {
              if (c == 0) {
                if (!line_data.fail()) {
                  gin_pts[b] = gin;
                } else {
                  Insist(0, "CSK data read failed!");
                }
              }
              // store the gamma out points for the first gamma in:
              if (b == 0) {
                if (!line_data.fail()) {
                  gout_pts[c] = gout;
                } else {
                  Insist(0, "CSK data read failed!");
                }
              }
            }

            for (size_t d = 0; d < n_xi; d++) // for each xi value
            {
              double csk_value;
              // get the data value and store it!
              line_data >> csk_value;
              if (!line_data.fail()) {
                raw_csk_data[a][k * gin_pts.size() + b][c][d] = csk_value;
              } else {
                Insist(0, "CSK data read failed!");
              }

            } // end xi loop
          }   // end gamma out loop
        }     // end gamma in loop
        // discard the blank line...
        std::string empty_line;
        getline(csk_data, empty_line);
      } // end loop over separate interpolation regions
    }   // end electron temperature loop

    Cdata->set_evalpts(etemp_pts, gin_pts, gout_pts, xi_pts);
    Cdata->set_csk_data(raw_csk_data);
  }
  // *************************************************************************
  // ************************** CLOSE THE FILE *******************************
  // *************************************************************************
  csk_data.close();

  std::cout << "CSK data read successfully!" << std::endl;

  return Cdata;
}

rtt_dsxx::SP<ComptonData> ComptonFile::read_lagrange_binary_csk_data() {
  // *************************************************************************
  // ************************ OPEN THE DATA FILE *****************************
  // *************************************************************************
  // open the file:
  csk_data.open(libfile, std::ios::binary);

  // *************************************************************************
  // ********************* READ SIZES/ALLOCATE ARRAYS ************************
  // *************************************************************************
  // first line of the file tells us how to allocate the raw data array
  // (which is the return value of the function:)
  if (!csk_data.is_open()) {
    std::ostringstream msg;
    msg << "ACK, we failed to open " << libfile << "! ";
    Insist(0, msg.str());
  }

  // *************************************************************************
  // ************************* GET DATA ARRAY SIZES **************************
  // *************************************************************************
  // declare variables for number of eval pts in gamma in, gamma out, xi, and
  // etemp:
  size_t n_ginbp, n_goutbp, n_etempbp;
  size_t n_gin, n_gout, n_xi, n_etemp;

  // make a vector to hold the raw size data:
  std::vector<char> size_data(7 * sizeof(size_t));

  // grab the first value (number of electron temp points)
  csk_data.read(&size_data[0], 7 * sizeof(size_t));

  // check for failure
  if (csk_data.fail()) {
    Insist(0, "Failed to read CSK data sizes!");
  }

  // cast raw character data to size_t
  std::memcpy(&n_etempbp, &size_data[0], sizeof(size_t));
  std::memcpy(&n_etemp, &size_data[sizeof(size_t)], sizeof(size_t));

  // cast raw character data to size_t
  std::memcpy(&n_ginbp, &size_data[2 * sizeof(size_t)], sizeof(size_t));
  std::memcpy(&n_gin, &size_data[3 * sizeof(size_t)], sizeof(size_t));

  // cast raw character data to size_t
  std::memcpy(&n_goutbp, &size_data[4 * sizeof(size_t)], sizeof(size_t));
  std::memcpy(&n_gout, &size_data[5 * sizeof(size_t)], sizeof(size_t));

  // cast raw character data to size_t
  std::memcpy(&n_xi, &size_data[6 * sizeof(size_t)], sizeof(size_t));

  // initialize data container for the "raw" csk values:
  SP_CompData Cdata(new ComptonData(n_etemp, n_gin, n_gout, n_xi));

  // *************************************************************************
  // *********************** GET THE BREAKPOINTS *****************************
  // *************************************************************************
  // get the next n_xi * sizeof(double) bytes from the file
  std::vector<double> etemp_breakpts(n_etempbp, 0.0);
  std::vector<char> etemp_data(n_etempbp * sizeof(double));
  csk_data.read(&etemp_data[0], n_etempbp * sizeof(double));
  // check for failure:
  if (csk_data.fail()) {
    Insist(0, "Failed to read electron temperature breakpoints!");
  }

  for (size_t m = 0; m < n_etempbp; m++) {
    double etempbp;
    std::memcpy(&etempbp, &etemp_data[m * sizeof(double)], sizeof(double));
    etemp_breakpts[m] = etempbp;
  }

  // get the next n_xi * sizeof(double) bytes from the file
  std::vector<double> gin_breakpts(n_ginbp, 0.0);
  std::vector<char> gin_data(n_ginbp * sizeof(double));
  csk_data.read(&gin_data[0], n_ginbp * sizeof(double));
  // check for failure:
  if (csk_data.fail()) {
    Insist(0, "Failed to read gin breakpoints!");
  }

  for (size_t m = 0; m < n_ginbp; m++) {
    double ginbp;
    std::memcpy(&ginbp, &gin_data[m * sizeof(double)], sizeof(double));
    gin_breakpts[m] = ginbp;
  }

  // get the next n_xi * sizeof(double) bytes from the file
  std::vector<double> gout_breakpts(n_goutbp, 0.0);
  std::vector<char> gout_data(n_goutbp * sizeof(double));
  csk_data.read(&gout_data[0], n_goutbp * sizeof(double));
  // check for failure:
  if (csk_data.fail()) {
    Insist(0, "Failed to read gout breakpoints!");
  }

  for (size_t m = 0; m < n_goutbp; m++) {
    double goutbp;
    std::memcpy(&goutbp, &gout_data[m * sizeof(double)], sizeof(double));
    gout_breakpts[m] = goutbp;
  }

  Cdata->set_breakpts(etemp_breakpts, gin_breakpts, gout_breakpts);

  // *************************************************************************
  // ************************ GET THE XI POINTS ******************************
  // *************************************************************************

  // get the next n_xi * sizeof(double) bytes from the file
  std::vector<double> xi_pts(n_xi, 0.0);
  std::vector<char> xi_data(n_xi * sizeof(double));
  csk_data.read(&xi_data[0], n_xi * sizeof(double));
  // check for failure:
  if (csk_data.fail()) {
    Insist(0, "Failed to read xi eval points!");
  }

  for (size_t m = 0; m < n_xi; m++) {
    double xi;
    std::memcpy(&xi, &xi_data[m * sizeof(double)], sizeof(double));
    xi_pts[m] = xi;
  }

  // *************************************************************************
  // ************************** READ THE DATA ********************************
  // *************************************************************************
  std::vector<double> etemp_pts(n_etemp, 0.0);
  std::vector<double> gin_pts(n_gin, 0.0);
  std::vector<double> gout_pts(n_gout, 0.0);

  std::vector<std::vector<std::vector<std::vector<double>>>> raw_csk_data(
      n_etemp, std::vector<std::vector<std::vector<double>>>(
                   3 * n_gin, std::vector<std::vector<double>>(
                                  n_gout, std::vector<double>(n_xi, 0.0))));

  for (size_t a = 0; a < n_etemp; a++) // for each electron temp expected
  {
    // get the electron temp for this block:
    double etemp;
    std::vector<char> etemp_val(sizeof(double));
    csk_data.read(&etemp_val[0], sizeof(double));
    memcpy(&etemp, &etemp_val[0], sizeof(double));
    etemp_pts[a] = etemp;

    // loop over the three regions in gin/gout, divided by boundary layers
    for (int k = 0; k < 3; k++) {
      // get the data for the block:
      std::vector<char> data_block(n_gin * n_gout * (n_xi + 2) *
                                   sizeof(double));
      csk_data.read(&data_block[0],
                    n_gin * n_gout * (n_xi + 2) * sizeof(double));

      for (size_t b = 0; b < n_gin; b++) // for each gamma in
      {
        for (size_t c = 0; c < n_gout; c++) // for each gamma out
        {
          // calculate the offset into the data array
          size_t offset = (b * n_gout + c) * (n_xi + 2) * sizeof(double);
          // get the first two values from the data block:
          double gin, gout;
          memcpy(&gin, &data_block[offset], sizeof(double));
          memcpy(&gout, &data_block[offset + sizeof(double)], sizeof(double));
          // store the gamma in points for the first electron temperature:
          if (a == 0) {
            if (c == 0) {
              gin_pts[b] = gin;
            }
            // store the gamma out points for the first gamma in:
            if (b == 0) {
              gout_pts[c] = gout;
            }
          }
          // update the data offset
          offset += 2 * sizeof(double);
          for (size_t d = 0; d < n_xi; d++) // for each xi value
          {
            double csk_value;
            // copy from the data array
            std::memcpy(&csk_value, &data_block[offset + d * sizeof(double)],
                        sizeof(double));
            raw_csk_data[a][k * gin_pts.size() + b][c][d] = csk_value;

          } // end xi loop
        }   // end gamma out loop
      }     // end gamma in loop
    }
  } // end electron temperature loop

  Cdata->set_evalpts(etemp_pts, gin_pts, gout_pts, xi_pts);
  Cdata->set_csk_data(raw_csk_data);

  // *************************************************************************
  // ************************** CLOSE THE FILE *******************************
  // *************************************************************************
  csk_data.close();

  std::cout << "CSK data read successfully!" << std::endl;
  return Cdata;
}
}
