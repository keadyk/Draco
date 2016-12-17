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

// prototype functions to check data
void check_mg_points(rtt_dsxx::ScalarUnitTest &, const std::vector<double> &,
                     const std::vector<double> &, const std::vector<double> &);

void check_mg_data(
    rtt_dsxx::ScalarUnitTest &,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &);

void tst_asciiread(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "======= Testing ascii CSK library read =======" << std::endl;
  std::cout << "==============================================" << std::endl;
  // ----------------------------------------------- //
  // Test the data file "mg_ascii"                  //
  // ----------------------------------------------- //
  std::string filename = "mg_ascii.compton";

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
  SP<ComptonData> spData = spCompton->read_mg_data();

  std::vector<std::vector<std::vector<std::vector<double>>>> tst_data =
      spData->get_csk_data();
  std::vector<double> gin_pts = spData->get_group_bds();
  std::vector<double> etemp_pts = spData->get_etemp_pts();
  std::vector<double> etemp_bpts = spData->get_etemp_breakpts();

  // check the points to be sure they match the file:
  check_mg_points(ut, gin_pts, etemp_pts, etemp_bpts);
  check_mg_data(ut, tst_data);
}

void check_mg_points(rtt_dsxx::ScalarUnitTest &ut,
                     const std::vector<double> &grp_bds,
                     const std::vector<double> &etemps,
                     const std::vector<double> &etemp_bps) {

  if (grp_bds.size() != 2) {
    ITFAILS;
  }
  if (etemps.size() != 7) {
    ITFAILS;
  }
  if (etemp_bps.size() != 2) {
    ITFAILS;
  }

  if (!soft_equiv(grp_bds[0], 3.91389000e-02, 1e-5))
    ITFAILS;
  if (!soft_equiv(grp_bds[1], 5.87084000e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(etemps[0], 1.76377944e-05, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemps[1], 8.95781632e-05, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemps[2], 2.05917691e-04, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemps[3], 3.46572429e-04, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemps[4], 4.87227167e-04, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemps[5], 6.03566703e-04, 1e-5))
    ITFAILS;
  if (!soft_equiv(etemps[6], 6.75507064e-04, 1e-5))
    ITFAILS;

  if (!soft_equiv(etemp_bps[0], 1.00000000e-06))
    ITFAILS;
  if (!soft_equiv(etemp_bps[1], 7.00000000e-04, 1e-5))
    ITFAILS;
}

void check_mg_data(
    rtt_dsxx::ScalarUnitTest &ut,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &data) {
  if (!soft_equiv(data[0][0][0][0], 4.47448782e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][0][1], 3.31200919e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][0][2], 4.53843561e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[0][0][0][3], 3.61765701e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(data[1][0][0][0], 4.47463414e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[1][0][0][1], 3.29459545e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[1][0][0][2], 4.53064161e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[1][0][0][3], 3.57345992e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(data[2][0][0][0], 4.47215388e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[2][0][0][1], 3.25531115e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[2][0][0][2], 4.51774000e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[2][0][0][3], 3.54424408e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(data[3][0][0][0], 4.46572469e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[3][0][0][1], 3.20974944e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[3][0][0][2], 4.50774880e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[3][0][0][3], 3.55958916e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(data[4][0][0][0], 4.45668383e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[4][0][0][1], 3.17337784e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[4][0][0][2], 4.50133379e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[4][0][0][3], 3.59663442e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(data[5][0][0][0], 4.44786593e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[5][0][0][1], 3.15048934e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[5][0][0][2], 4.49718508e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[5][0][0][3], 3.63432689e-02, 1e-5))
    ITFAILS;

  if (!soft_equiv(data[6][0][0][0], 4.44197781e+00, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[6][0][0][1], 3.13926441e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[6][0][0][2], 4.49480847e-01, 1e-5))
    ITFAILS;
  if (!soft_equiv(data[6][0][0][3], 3.65912346e-02, 1e-5))
    ITFAILS;
}

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_asciiread(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
