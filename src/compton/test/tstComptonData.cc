//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/tstComptonData.cc
 * \author Kendra P. Keady
 * \date   Mon Oct 17 2016
 * \brief  Compton data class test.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */

#include "ComptonData.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/SP.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>

using rtt_compton::ComptonData;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

void tst_accessors(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "====== Testing ComptonData construction ======" << std::endl;
  std::cout << "==============================================" << std::endl;

  size_t n_etempbp = 5;

  size_t n_etemp = 2;
  size_t n_grp = 33;
  size_t n_leg = 41;

  double etemp_bpval = 8888.0;

  double etemp_val = 5555.0;
  double grp_val = 4444.0;

  // create data arrays:
  std::vector<double> etemp_bpts(n_etempbp, etemp_bpval);

  std::vector<double> etemp_pts(n_etemp, etemp_val);
  std::vector<double> grp_bds(n_grp, grp_val);

  ComptonData Cdata(n_etemp, n_grp, n_leg);

  // set data arrays in the Cdata object:
  Cdata.set_evalpts(etemp_pts, grp_bds);
  Cdata.set_breakpts(etemp_bpts);

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      n_etemp, std::vector<std::vector<std::vector<double>>>(
                   n_grp, std::vector<std::vector<double>>(
                              n_grp, std::vector<double>(n_leg, 0.0))));
  fake_data[0][0][0][0] = 1.2345;
  fake_data[n_etemp - 1][n_grp - 1][n_grp - 1][n_leg - 1] = 6.789;
  Cdata.set_csk_data(fake_data);

  // test data sizes:
  if (Cdata.get_n_etemp_breakpts() != n_etempbp)
    ITFAILS;

  if (Cdata.get_n_etemp_pts() != n_etemp)
    ITFAILS;

  if (Cdata.get_n_leg() != n_leg)
    ITFAILS;

  if (Cdata.get_csk_data()[0][0][0][0] != 1.2345)
    ITFAILS;
  if (Cdata.get_csk_data()[n_etemp - 1][n_grp - 1][n_grp - 1][n_leg - 1] !=
      6.789)
    ITFAILS;

  // test all of the non-const get routines:
  {
    std::vector<double> etemp_bpcomp = Cdata.get_etemp_breakpts();
    std::vector<double> etemp_comp = Cdata.get_etemp_pts();
    std::vector<double> grp_bds_comp = Cdata.get_group_bds();

    for (size_t a = 0; a < etemp_bpcomp.size(); a++) {
      if (etemp_bpcomp[a] != etemp_bpval)
        ITFAILS;
    }

    for (size_t a = 0; a < etemp_comp.size(); a++) {
      if (etemp_comp[a] != etemp_val)
        ITFAILS;
    }
    for (size_t a = 0; a < grp_bds_comp.size(); a++) {
      if (grp_bds_comp[a] != grp_val)
        ITFAILS;
    }
  }

  if (ut.numFails == 0) {
    PASSMSG("All non-const ComptonData accessors functioned correctly.");
  } else {
    FAILMSG("Some non-const ComptonData accessors functioned incorrectly!");
  }

  const ComptonData *cCdata(&Cdata);

  // test data sizes:
  if (cCdata->get_n_etemp_breakpts() != n_etempbp)
    ITFAILS;

  if (cCdata->get_n_etemp_pts() != n_etemp)
    ITFAILS;

  if (cCdata->get_n_leg() != n_leg)
    ITFAILS;

  if (cCdata->get_csk_data()[0][0][0][0] != 1.2345)
    ITFAILS;
  if (cCdata->get_csk_data()[n_etemp - 1][n_grp - 1][n_grp - 1][n_leg - 1] !=
      6.789)
    ITFAILS;

  // test all of the const get routines:
  {
    std::vector<double> etemp_bpcomp = cCdata->get_etemp_breakpts();

    std::vector<double> etemp_comp = cCdata->get_etemp_pts();
    std::vector<double> grp_bds_comp = cCdata->get_group_bds();

    for (size_t a = 0; a < etemp_bpcomp.size(); a++) {
      if (etemp_bpcomp[a] != etemp_bpval)
        ITFAILS;
    }
  }

  if (ut.numFails == 0) {
    PASSMSG("All const ComptonData accessors functioned correctly.");
  } else {
    FAILMSG("Some const ComptonData accessors functioned incorrectly!");
  }

  if (ut.numFails == 0) {
    std::cout << "ComptonData accessor unit test PASSED." << std::endl;
  } else {
    std::cout << "ComptonData accessor unit test FAILED." << std::endl;
  }
}

void tst_mg_data(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "=== Testing multigrp ComptonData container ===" << std::endl;
  std::cout << "==============================================" << std::endl;
  // Since ComptonData is a very simple data container with setters and getters,
  // all we can really do is plug in weird numbers and make sure we get them
  // back...

  // make up some fake data...
  std::vector<double> etemp(2, 0.0);
  etemp[0] = 0.056;
  etemp[1] = 1.00894;
  std::vector<double> grp_bds(4, 0.0);
  grp_bds[0] = 0.0005;
  grp_bds[1] = 0.54321;
  grp_bds[2] = 0.8974;
  grp_bds[2] = 0.99;
  size_t nleg = 5;

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      2,
      std::vector<std::vector<std::vector<double>>>(
          3, std::vector<std::vector<double>>(3, std::vector<double>(5, 0.0))));

  fake_data[0][0][0][0] = 1.345e-15;
  fake_data[0][0][1][1] = 8.888e-46;
  fake_data[0][1][0][2] = 2.687e-01;
  fake_data[0][1][1][3] = 8.888e-04;
  fake_data[1][1][0][4] = 1.000;
  fake_data[1][1][1][0] = 3.141e-29;
  fake_data[1][2][0][1] = 9.999e-99;
  fake_data[1][2][1][2] = 4.646e-02;

  // try creating a compton data object:
  SP<ComptonData> Cdata;
  try {
    Cdata.reset(new ComptonData(etemp.size(), (grp_bds.size() - 1), nleg));
    PASSMSG("ComptonData object successfully created!");
  } catch (rtt_dsxx::assertion &as) {
    FAILMSG(as.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonData object";
    FAILMSG(message.str());
    return;
  }

  // next, set all of the standard data vectors:
  Cdata->set_evalpts(etemp, grp_bds);
  Cdata->set_csk_data(fake_data);

  std::vector<double> etemp_check = Cdata->get_etemp_pts();
  std::vector<double> grp_bds_check = Cdata->get_group_bds();
  std::vector<std::vector<std::vector<std::vector<double>>>> data_check =
      Cdata->get_csk_data();

  // check the sizes
  if (etemp.size() != etemp_check.size())
    ITFAILS;
  if (grp_bds.size() != grp_bds_check.size())
    ITFAILS;

  if (fake_data.size() != data_check.size())
    ITFAILS;
  if (fake_data[0].size() != data_check[0].size())
    ITFAILS;
  if (fake_data[0][0].size() != data_check[0][0].size())
    ITFAILS;
  if (fake_data[0][0][0].size() != data_check[0][0][0].size())
    ITFAILS;

  // Now, scroll through and check all the get() routines against the values
  // we passed to set():
  for (size_t a = 0; a < etemp.size(); a++) {
    if (!soft_equiv(etemp[a], etemp_check[a]))
      ITFAILS;
    for (size_t b = 0; b < (grp_bds.size() - 1); b++) {
      for (size_t c = 0; c < (grp_bds.size() - 1); c++) {
        for (size_t d = 0; d < nleg; d++) {
          if (!soft_equiv(data_check[a][b][c][d], fake_data[a][b][c][d]))
            ITFAILS;
        }
      }
    }
  }

  // check grp_bds, gout, and xi points:
  for (size_t b = 0; b < grp_bds.size(); b++) {
    if (!soft_equiv(grp_bds[b], grp_bds_check[b]))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("All ComptonData CSK values returned correctly.");
  } else {
    FAILMSG("Some ComptonData CSK values incorrect!");
  }

  if (ut.numFails == 0) {
    std::cout << "Standard ComptonData container unit test PASSED."
              << std::endl;
  } else {
    std::cout << "Standard ComptonData container unit test FAILED."
              << std::endl;
  }
}

void tst_data_slices(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "== Testing ComptonData data subset routines ==" << std::endl;
  std::cout << "==============================================" << std::endl;
  // make up some fake data...
  // (we tested get() and set() for the eval points in the standard test,
  // so here we only put data in the breakpoint vectors)
  std::vector<double> etemp_bp(3, 0.0);
  etemp_bp[0] = 0.3333;
  etemp_bp[1] = 0.6666;
  etemp_bp[2] = 1.9999;
  std::vector<double> grp_bds(5, 0.0);
  grp_bds[0] = 0.25;
  grp_bds[1] = 0.35;
  grp_bds[2] = 0.65;
  grp_bds[3] = 0.75;
  grp_bds[4] = 0.88;
  std::vector<double> etemp(2, 0.0);
  etemp[0] = 0.4444;
  etemp[1] = 0.7777;
  size_t nleg = 1;

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      2,
      std::vector<std::vector<std::vector<double>>>(
          4, std::vector<std::vector<double>>(4, std::vector<double>(1, 0.0))));

  // fill in the data with simple increasing doubles
  for (size_t a = 0; a < fake_data.size(); a++) {
    for (size_t b = 0; b < fake_data[a].size(); b++) {
      for (size_t c = 0; c < fake_data[a][b].size(); c++) {
        for (size_t d = 0; d < fake_data[a][b][c].size(); d++) {
          fake_data[a][b][c][d] =
              static_cast<double>(a * (12 * 4) + b * 4 + c) + 0.5;
        }
      }
    }
  }

  // try creating a compton data object:
  SP<ComptonData> Cdata;
  try {
    Cdata.reset(new ComptonData(etemp.size(), (grp_bds.size() - 1), nleg));
    PASSMSG("ComptonData object successfully created!");
  } catch (rtt_dsxx::assertion &as) {
    FAILMSG(as.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonData object";
    FAILMSG(message.str());
    return;
  }

  // next, set all of the standard and breakpoint vectors in the ComptonData
  // container:
  Cdata->set_breakpts(etemp_bp);
  Cdata->set_evalpts(etemp, grp_bds);
  Cdata->set_csk_data(fake_data);

  // get the slice of data corresponding to electron temperature index 0:
  std::vector<std::vector<std::vector<std::vector<double>>>> slice0 =
      Cdata->get_csk_data(0);

  // get the slice of data corresponding to electron temperature index 1:
  std::vector<std::vector<std::vector<std::vector<double>>>> slice1 =
      Cdata->get_csk_data(1);

  // check the returned data against the origrp_bdsal values
  for (size_t a = 0; a < slice1.size(); a++) {
    for (size_t b = 0; b < slice1[a].size(); b++) {
      for (size_t c = 0; c < slice1[a][b].size(); c++) {
        for (size_t d = 0; d < slice1[a][b][c].size(); d++) {
          if (!soft_equiv(slice0[a][b][c][d], fake_data[0][b][c][d])) {
            ITFAILS;
          }
          if (!soft_equiv(slice1[a][b][c][d], fake_data[1][b][c][d])) {
            ITFAILS;
          }
        }
      }
    }
  }
  if (ut.numFails == 0) {
    PASSMSG("All ComptonData large data slices returned correctly.");
  } else {
    FAILMSG("Some ComptonData large data slice values incorrect!");
  }
}

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_accessors(ut);
    tst_mg_data(ut);
    tst_data_slices(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
