//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/tstComptonFile.cc
 * \author Kendra P. Keady
 * \date   Mon Oct 17 2016
 * \brief  Compton test.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
#include "ComptonData.hh"
#include "ComptonFile.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/SP.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>

using rtt_compton::ComptonFile;
using rtt_compton::ComptonData;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

// prototype functions to check data:
void check_std_points(rtt_dsxx::ScalarUnitTest &, const std::vector<double> &,
                      const std::vector<double> &, const std::vector<double> &,
                      const std::vector<double> &);
void check_lag_points(rtt_dsxx::ScalarUnitTest &, const std::vector<double> &,
                      const std::vector<double> &, const std::vector<double> &,
                      const std::vector<double> &, const std::vector<double> &,
                      const std::vector<double> &, const std::vector<double> &);

void check_std_data(
    rtt_dsxx::ScalarUnitTest &,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &);

void check_lag_data(
    rtt_dsxx::ScalarUnitTest &,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &);

void tst_asciiread(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "======= Testing ascii CSK library read =======" << std::endl;
  std::cout << "==============================================" << std::endl;
  // ----------------------------------------------- //
  // Test the data file "csk_ascii"                  //
  // ----------------------------------------------- //
  std::string filename = "csk_ascii.compton";

  // Create a smart pointer to a ComptonFile object
  SP<ComptonFile> spCompton;

  // Try to instantiate the object.
  try {
    spCompton.reset(new ComptonFile(filename, false));
  } catch (rtt_dsxx::assertion const &excpt) {
    FAILMSG(excpt.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonFile object";
    FAILMSG(message.str());
    return;
  }

  // Success! Carry on, carry on.
  PASSMSG("Successfully created ComptonFile object.");

  // Next, read in the data
  SP<ComptonData> spData = spCompton->read_csk_data();

  std::vector<std::vector<std::vector<std::vector<double>>>> tst_data =
      spData->get_csk_data();
  std::vector<double> gin_pts = spData->get_gin_pts();
  std::vector<double> gout_pts = spData->get_gout_pts();
  std::vector<double> xi_pts = spData->get_xi_pts();
  std::vector<double> etemp_pts = spData->get_etemp_pts();

  // check the points to be sure they match the file:
  check_std_points(ut, gin_pts, gout_pts, xi_pts, etemp_pts);
  check_std_data(ut, tst_data);
}

void tst_binread(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "===============================================" << std::endl;
  std::cout << "======= Testing binary CSK library read =======" << std::endl;
  std::cout << "===============================================" << std::endl;
  // ----------------------------------------------- //
  // Test the data file "csk_bin"                    //
  // ----------------------------------------------- //
  std::string filename = "csk_bin.compton";

  // Create a smart pointer to a ComptonFile object
  SP<ComptonFile> spCompton;

  // Try to instantiate the object.
  try {
    spCompton.reset(new ComptonFile(filename));
  } catch (rtt_dsxx::assertion const &excpt) {
    FAILMSG(excpt.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonFile object";
    FAILMSG(message.str());
    return;
  }

  // Success! Carry on, carry on.
  PASSMSG("Successfully created ComptonFile object.");

  // Next, read in the data
  SP<ComptonData> spData = spCompton->read_csk_data();

  std::vector<std::vector<std::vector<std::vector<double>>>> tst_data =
      spData->get_csk_data();
  std::vector<double> gin_pts = spData->get_gin_pts();
  std::vector<double> gout_pts = spData->get_gout_pts();
  std::vector<double> xi_pts = spData->get_xi_pts();
  std::vector<double> etemp_pts = spData->get_etemp_pts();

  // check the points to be sure they match the file:
  check_std_points(ut, gin_pts, gout_pts, xi_pts, etemp_pts);
  check_std_data(ut, tst_data);
}

void tst_lagrange_asciiread(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "===============================================" << std::endl;
  std::cout << "=== Testing ascii CSK Lagrange library read ===" << std::endl;
  std::cout << "===============================================" << std::endl;
  // ----------------------------------------------- //
  // Test the data file "lagrange_csk_ascii"         //
  // ----------------------------------------------- //
  std::string filename = "lagrange_csk_ascii.compton";

  // Create a smart pointer to a ComptonFile object
  SP<ComptonFile> spCompton;

  // Try to instantiate the object.
  try {
    spCompton.reset(new ComptonFile(filename, false));
  } catch (rtt_dsxx::assertion const &excpt) {
    FAILMSG(excpt.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonFile object";
    FAILMSG(message.str());
    return;
  }

  // Success! Carry on, carry on.
  PASSMSG("Successfully created ComptonFile object.");

  // Next, read in the data
  SP<ComptonData> spData = spCompton->read_lagrange_csk_data();

  std::vector<std::vector<std::vector<std::vector<double>>>> tst_data =
      spData->get_csk_data();
  std::vector<double> gin_pts = spData->get_gin_pts();
  std::vector<double> gin_breakpts = spData->get_gin_breakpts();
  std::vector<double> gout_pts = spData->get_gout_pts();
  std::vector<double> gout_breakpts = spData->get_gout_breakpts();
  std::vector<double> xi_pts = spData->get_xi_pts();
  std::vector<double> etemp_pts = spData->get_etemp_pts();
  std::vector<double> etemp_breakpts = spData->get_etemp_breakpts();

  check_lag_points(ut, gin_breakpts, gin_pts, gout_breakpts, gout_pts, xi_pts,
                   etemp_breakpts, etemp_pts);
  check_lag_data(ut, tst_data);
}

void tst_lagrange_binread(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "================================================" << std::endl;
  std::cout << "=== Testing binary CSK Lagrange library read ===" << std::endl;
  std::cout << "================================================" << std::endl;
  // ----------------------------------------------- //
  // Test the data file "lagrange_csk_binary"        //
  // ----------------------------------------------- //
  std::string filename = "lagrange_csk_binary.compton";

  // Create a smart pointer to a ComptonFile object
  SP<ComptonFile> spCompton;

  // Try to instantiate the object.
  try {
    spCompton.reset(new ComptonFile(filename));
  } catch (rtt_dsxx::assertion const &excpt) {
    FAILMSG(excpt.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonFile object";
    FAILMSG(message.str());
    return;
  }

  // Success! Carry on, carry on.
  PASSMSG("Successfully created ComptonFile object.");

  // Next, read in the data
  SP<ComptonData> spData = spCompton->read_lagrange_csk_data();

  std::vector<std::vector<std::vector<std::vector<double>>>> tst_data =
      spData->get_csk_data();

  std::vector<double> gin_pts = spData->get_gin_pts();
  std::vector<double> gin_breakpts = spData->get_gin_breakpts();
  std::vector<double> gout_pts = spData->get_gout_pts();
  std::vector<double> gout_breakpts = spData->get_gout_breakpts();
  std::vector<double> xi_pts = spData->get_xi_pts();
  std::vector<double> etemp_pts = spData->get_etemp_pts();
  std::vector<double> etemp_breakpts = spData->get_etemp_breakpts();

  check_lag_points(ut, gin_breakpts, gin_pts, gout_breakpts, gout_pts, xi_pts,
                   etemp_breakpts, etemp_pts);
  check_lag_data(ut, tst_data);
}

// point check for standard library:
void check_std_points(rtt_dsxx::ScalarUnitTest &ut,
                      const std::vector<double> &gin_pts,
                      const std::vector<double> &gout_pts,
                      const std::vector<double> &xi_pts,
                      const std::vector<double> &etemp_pts) {
  double gin_ref, etemp_ref;
  // check the uniformly-discretized variables:
  for (size_t a = 0; a < 20; a++) {
    gin_ref = 0.05 + a * 0.1;
    if (!soft_equiv(gin_pts[a], gin_ref))
      ITFAILS;
  }
  // finely-discretized region
  for (size_t b = 0; b < 5; b++) {
    etemp_ref = 0.025 + b * 0.05;
    if (!soft_equiv(etemp_pts[b], etemp_ref))
      ITFAILS;
  }
  // coarsely-discretized region
  for (size_t c = 0; c < 5; c++) {
    etemp_ref = 0.325 + c * 0.15;
    if (!soft_equiv(etemp_pts[5 + c], etemp_ref))
      ITFAILS;
  }

  // check the arbitrarily-discretized data
  if (!soft_equiv(gout_pts[0], 0.01))
    ITFAILS;
  if (!soft_equiv(gout_pts[1], 0.035))
    ITFAILS;
  if (!soft_equiv(gout_pts[2], 0.28))
    ITFAILS;
  if (!soft_equiv(gout_pts[3], 0.5))
    ITFAILS;
  if (!soft_equiv(gout_pts[4], 0.99))
    ITFAILS;

  if (!soft_equiv(xi_pts[0], -0.95))
    ITFAILS;
  if (!soft_equiv(xi_pts[1], 0.0))
    ITFAILS;
  if (!soft_equiv(xi_pts[2], 0.25))
    ITFAILS;
  if (!soft_equiv(xi_pts[3], 0.88))
    ITFAILS;
  if (!soft_equiv(xi_pts[4], 0.95))
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully read evaluation points from test library!");
  } else {
    FAILMSG("Did not read correct evaluation points from test library!");
  }
}

// point check for standard library:
void check_lag_points(rtt_dsxx::ScalarUnitTest &ut,
                      const std::vector<double> &gin_breakpts,
                      const std::vector<double> &gin_pts,
                      const std::vector<double> &gout_breakpts,
                      const std::vector<double> &gout_pts,
                      const std::vector<double> &xi_pts,
                      const std::vector<double> &etemp_breakpts,
                      const std::vector<double> &etemp_pts) {
  // check the data sizes:
  if (gin_breakpts.size() != 3)
    ITFAILS;
  if (gin_pts.size() != 4)
    ITFAILS;
  if (gout_breakpts.size() != 3)
    ITFAILS;
  if (gout_pts.size() != 4)
    ITFAILS;
  if (etemp_breakpts.size() != 5)
    ITFAILS;
  if (etemp_pts.size() != 8)
    ITFAILS;
  if (xi_pts.size() != 5)
    ITFAILS;

  // check the global breakpoints
  if (!soft_equiv(gin_breakpts[0], 0.0))
    ITFAILS;
  if (!soft_equiv(gin_breakpts[1], 0.5))
    ITFAILS;
  if (!soft_equiv(gin_breakpts[2], 1.0))
    ITFAILS;
  if (!soft_equiv(gout_breakpts[0], 0.0))
    ITFAILS;
  if (!soft_equiv(gout_breakpts[1], 0.5))
    ITFAILS;
  if (!soft_equiv(gout_breakpts[2], 1.0))
    ITFAILS;
  if (!soft_equiv(etemp_breakpts[0], 0.0))
    ITFAILS;
  if (!soft_equiv(etemp_breakpts[1], 0.0223607, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemp_breakpts[2], 0.5))
    ITFAILS;
  if (!soft_equiv(etemp_breakpts[3], 0.977639, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemp_breakpts[4], 1.0))
    ITFAILS;

  // check the local evaluation points:
  if (!soft_equiv(gin_pts[0], 0.105662, 1e-5))
    ITFAILS;
  if (!soft_equiv(gin_pts[1], 0.394338, 1e-5))
    ITFAILS;
  if (!soft_equiv(gin_pts[2], 0.605662, 1e-5))
    ITFAILS;
  if (!soft_equiv(gin_pts[3], 0.894338, 1e-5))
    ITFAILS;

  if (!soft_equiv(gout_pts[0], 0.105662, 1e-5))
    ITFAILS;
  if (!soft_equiv(gout_pts[1], 0.394338, 1e-5))
    ITFAILS;
  if (!soft_equiv(gout_pts[2], 0.605662, 1e-5))
    ITFAILS;
  if (!soft_equiv(gout_pts[3], 0.894338, 1e-5))
    ITFAILS;

  // check the uniformly-discretized angular data:
  if (!soft_equiv(xi_pts[0], -0.8))
    ITFAILS;
  if (!soft_equiv(xi_pts[1], -0.4))
    ITFAILS;
  if (!soft_equiv(xi_pts[2], 0.0))
    ITFAILS;
  if (!soft_equiv(xi_pts[3], 0.4))
    ITFAILS;
  if (!soft_equiv(xi_pts[4], 0.8))
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully read evaluation points from test library!");
  } else {
    FAILMSG("Did not read correct evaluation points from test library!");
  }
}

//data check for standard library:
void check_std_data(
    rtt_dsxx::ScalarUnitTest &ut,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &data) {
  // check selected data points from each regime:
  if (data.size() != 10)
    ITFAILS;
  if (data[0].size() != 20)
    ITFAILS;
  if (data[0][0].size() != 5)
    ITFAILS;
  if (data[0][0][0].size() != 5)
    ITFAILS;
  // SCALE small values and compare to one:
  if (!soft_equiv(data[0][0][0][0], 5.55794e-05, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][0][4] / 2.02745e-81, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][4][0] / 3.07426e-31, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][4][4] / 2.9701e-222, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][19][0][0] / 1.07261e-92, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][19][0][4], 0.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][19][4][0], 2.0226e-05, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][19][4][4] / 3.45633e-17, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][0][0][0], 0.755117, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][0][0][4], 0.00901667, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][0][4][0], 0.862336, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][0][4][4], 6.41622e-06, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][19][0][0], 3.62824e-06, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][19][0][4] / 2.31686e-23, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][19][4][0], 0.124555, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[9][19][4][4], 0.0684458, 1e-5))
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully read CSK data points from test library!");
  } else {
    FAILMSG("Did not read correct CSK data points from test library!");
  }
}

//data check for lagrange library:
void check_lag_data(
    rtt_dsxx::ScalarUnitTest &ut,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &data) {
  // check selected data points from each regime:
  if (data.size() != 8)
    ITFAILS;
  if (data[0].size() != 3 * 4)
    ITFAILS;
  if (data[0][0].size() != 4)
    ITFAILS;
  if (data[0][0][0].size() != 5)
    ITFAILS;
  // SCALE small values and compare to one:
  if (!soft_equiv(data[0][0][0][0] / 1.29629e-44, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][0][4] / 2.59185e-48, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][3][0] / 1.85086e-57, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][3][4] / 5.30928e-141, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][11][0][0] / 1.58091e-51, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][11][0][4] / 4.58293e-171, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][11][3][0] / 3.51495e-147, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][11][3][4] / 1.01014e-142, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][0][0][0] / 2.20299e-46, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][0][0][4] / 7.52173e-46, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][0][3][0] / 1.03222e-46, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][0][3][4] / 7.49007e-47, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][11][0][0] / 1.20859e-47, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][11][0][4] / 4.03163e-48, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][11][3][0] / 3.09407e-47, 1.0, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[7][11][3][4] / 3.53662e-47, 1.0, 1e-5))
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully read CSK data points from test library!");
  } else {
    FAILMSG("Did not read correct CSK data points from test library!");
  }
}

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_asciiread(ut);
    tst_binread(ut);
    tst_lagrange_asciiread(ut);
    tst_lagrange_binread(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
