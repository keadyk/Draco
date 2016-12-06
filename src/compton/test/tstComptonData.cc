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
  size_t n_ginbp = 3;
  size_t n_goutbp = 10;

  size_t n_etemp = 2;
  size_t n_gin = 33;
  size_t n_gout = 33;
  size_t n_xi = 41;

  double etemp_bpval = 8888.0;
  double gin_bpval = 7777.0;
  double gout_bpval = 6666.0;

  double etemp_val = 5555.0;
  double gin_val = 4444.0;
  double gout_val = 3333.0;
  double xi_val = -0.2222;
  bool lagrange = true;
  bool legendre = false;

  // create data arrays:
  std::vector<double> etemp_bpts(n_etempbp, etemp_bpval);
  std::vector<double> gin_bpts(n_ginbp, gin_bpval);
  std::vector<double> gout_bpts(n_goutbp, gout_bpval);

  std::vector<double> etemp_pts(n_etemp, etemp_val);
  std::vector<double> gin_pts(n_gin, gin_val);
  std::vector<double> gout_pts(n_gout, gout_val);
  std::vector<double> xi_pts(n_xi, xi_val);

  ComptonData Cdata(n_etemp, n_gin, n_gout, n_xi, lagrange, legendre);

  // test bool accessors
  if (Cdata.is_legendre())
    ITFAILS;
  if (!Cdata.is_lagrange())
    ITFAILS;

  // set data arrays in the Cdata object:
  Cdata.set_evalpts(etemp_pts, gin_pts, gout_pts, xi_pts);
  Cdata.set_breakpts(etemp_bpts, gin_bpts, gout_bpts);

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      n_etemp, std::vector<std::vector<std::vector<double>>>(
                   3 * n_gin, std::vector<std::vector<double>>(
                                  n_gout, std::vector<double>(n_xi, 0.0))));
  fake_data[0][0][0][0] = 1.2345;
  fake_data[n_etemp - 1][3 * n_gin - 1][n_gout - 1][n_xi - 1] = 6.789;
  Cdata.set_csk_data(fake_data);

  // test data sizes:
  if (Cdata.get_n_etemp_breakpts() != n_etempbp)
    ITFAILS;
  if (Cdata.get_n_gin_breakpts() != n_ginbp)
    ITFAILS;
  if (Cdata.get_n_gout_breakpts() != n_goutbp)
    ITFAILS;

  if (Cdata.get_n_etemp_pts() != n_etemp)
    ITFAILS;
  if (Cdata.get_n_gin_pts() != n_gin)
    ITFAILS;
  if (Cdata.get_n_gout_pts() != n_gout)
    ITFAILS;
  if (Cdata.get_n_xi_pts() != n_xi)
    ITFAILS;

  if (Cdata.get_csk_data()[0][0][0][0] != 1.2345)
    ITFAILS;
  if (Cdata.get_csk_data()[n_etemp - 1][3 * n_gin - 1][n_gout - 1][n_xi - 1] !=
      6.789)
    ITFAILS;

  // test all of the non-const get routines:
  {
    std::vector<double> etemp_bpcomp = Cdata.get_etemp_breakpts();
    std::vector<double> gin_bpcomp = Cdata.get_gin_breakpts();
    std::vector<double> gout_bpcomp = Cdata.get_gout_breakpts();

    std::vector<double> etemp_comp = Cdata.get_etemp_pts();
    std::vector<double> gin_comp = Cdata.get_gin_pts();
    std::vector<double> gout_comp = Cdata.get_gout_pts();
    std::vector<double> xi_comp = Cdata.get_xi_pts();

    for (size_t a = 0; a < etemp_bpcomp.size(); a++) {
      if (etemp_bpcomp[a] != etemp_bpval)
        ITFAILS;
    }
    for (size_t a = 0; a < gin_bpcomp.size(); a++) {
      if (gin_bpcomp[a] != gin_bpval)
        ITFAILS;
    }
    for (size_t a = 0; a < gout_bpcomp.size(); a++) {
      if (gout_bpcomp[a] != gout_bpval)
        ITFAILS;
    }
    for (size_t a = 0; a < etemp_comp.size(); a++) {
      if (etemp_comp[a] != etemp_val)
        ITFAILS;
    }
    for (size_t a = 0; a < gin_comp.size(); a++) {
      if (gin_comp[a] != gin_val)
        ITFAILS;
    }
    for (size_t a = 0; a < gout_comp.size(); a++) {
      if (gout_comp[a] != gout_val)
        ITFAILS;
    }
    for (size_t a = 0; a < xi_comp.size(); a++) {
      if (xi_comp[a] != xi_val)
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
  if (cCdata->get_n_gin_breakpts() != n_ginbp)
    ITFAILS;
  if (cCdata->get_n_gout_breakpts() != n_goutbp)
    ITFAILS;

  if (cCdata->get_n_etemp_pts() != n_etemp)
    ITFAILS;
  if (cCdata->get_n_gin_pts() != n_gin)
    ITFAILS;
  if (cCdata->get_n_gout_pts() != n_gout)
    ITFAILS;
  if (cCdata->get_n_xi_pts() != n_xi)
    ITFAILS;

  if (cCdata->get_csk_data()[0][0][0][0] != 1.2345)
    ITFAILS;
  if (cCdata->get_csk_data()[n_etemp - 1][3 * n_gin - 1][n_gout - 1]
                            [n_xi - 1] != 6.789)
    ITFAILS;

  // test all of the const get routines:
  {
    std::vector<double> etemp_bpcomp = cCdata->get_etemp_breakpts();
    std::vector<double> gin_bpcomp = cCdata->get_gin_breakpts();
    std::vector<double> gout_bpcomp = cCdata->get_gout_breakpts();

    std::vector<double> etemp_comp = cCdata->get_etemp_pts();
    std::vector<double> gin_comp = cCdata->get_gin_pts();
    std::vector<double> gout_comp = cCdata->get_gout_pts();
    std::vector<double> xi_comp = cCdata->get_xi_pts();

    for (size_t a = 0; a < etemp_bpcomp.size(); a++) {
      if (etemp_bpcomp[a] != etemp_bpval)
        ITFAILS;
    }
    for (size_t a = 0; a < gin_bpcomp.size(); a++) {
      if (gin_bpcomp[a] != gin_bpval)
        ITFAILS;
    }
    for (size_t a = 0; a < gout_bpcomp.size(); a++) {
      if (gout_bpcomp[a] != gout_bpval)
        ITFAILS;
    }
    for (size_t a = 0; a < etemp_comp.size(); a++) {
      if (etemp_comp[a] != etemp_val)
        ITFAILS;
    }
    for (size_t a = 0; a < gin_comp.size(); a++) {
      if (gin_comp[a] != gin_val)
        ITFAILS;
    }
    for (size_t a = 0; a < gout_comp.size(); a++) {
      if (gout_comp[a] != gout_val)
        ITFAILS;
    }
    for (size_t a = 0; a < xi_comp.size(); a++) {
      if (xi_comp[a] != xi_val)
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

void tst_std_data(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "=== Testing standard ComptonData container ===" << std::endl;
  std::cout << "==============================================" << std::endl;
  // Since ComptonData is a very simple data container with setters and getters,
  // all we can really do is plug in weird numbers and make sure we get them
  // back...

  // make up some fake data...
  std::vector<double> etemp(2, 0.0);
  etemp[0] = 0.056;
  etemp[1] = 1.00894;
  std::vector<double> gin(3, 0.0);
  gin[0] = 0.0005;
  gin[1] = 0.54321;
  gin[2] = 0.8974;
  std::vector<double> gout(2, 0.0);
  gout[0] = 0.1;
  gout[1] = 0.7;
  std::vector<double> xi(5, 0.0);
  xi[0] = -0.99;
  xi[1] = -0.45;
  xi[2] = 0.0;
  xi[3] = 0.45;
  xi[4] = 0.99;

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      2,
      std::vector<std::vector<std::vector<double>>>(
          3, std::vector<std::vector<double>>(2, std::vector<double>(5, 0.0))));

  fake_data[0][0][0][0] = 1.345e-15;
  fake_data[0][0][1][1] = 8.888e-46;
  fake_data[0][1][0][2] = 2.687e-01;
  fake_data[0][1][1][3] = 8.888e-04;
  fake_data[1][1][0][4] = 1.000;
  fake_data[1][1][1][0] = 3.141e-29;
  fake_data[1][2][0][1] = 9.999e-99;
  fake_data[1][2][1][2] = 4.646e-02;

  // set interpolation type bools:
  bool lagrange = false;
  bool legendre = true;

  // try creating a compton data object:
  SP<ComptonData> Cdata;
  try {
    Cdata.reset(new ComptonData(etemp.size(), gin.size(), gout.size(),
                                xi.size(), lagrange, legendre));
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
  Cdata->set_evalpts(etemp, gin, gout, xi);
  Cdata->set_csk_data(fake_data);

  std::vector<double> etemp_check = Cdata->get_etemp_pts();
  std::vector<double> gin_check = Cdata->get_gin_pts();
  std::vector<double> gout_check = Cdata->get_gout_pts();
  std::vector<double> xi_check = Cdata->get_xi_pts();
  std::vector<std::vector<std::vector<std::vector<double>>>> data_check =
      Cdata->get_csk_data();

  // check the sizes
  if (etemp.size() != etemp_check.size())
    ITFAILS;
  if (gin.size() != gin_check.size())
    ITFAILS;
  if (gout.size() != gout_check.size())
    ITFAILS;
  if (xi.size() != xi_check.size())
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
    for (size_t b = 0; b < gin.size(); b++) {
      for (size_t c = 0; c < gout.size(); c++) {
        for (size_t d = 0; d < xi.size(); d++) {
          if (!soft_equiv(data_check[a][b][c][d], fake_data[a][b][c][d]))
            ITFAILS;
        }
      }
    }
  }

  // check gin, gout, and xi points:
  for (size_t b = 0; b < gin.size(); b++) {
    if (!soft_equiv(gin[b], gin_check[b]))
      ITFAILS;
  }
  for (size_t c = 0; c < gout.size(); c++) {
    if (!soft_equiv(gout[c], gout_check[c]))
      ITFAILS;
  }
  for (size_t d = 0; d < xi.size(); d++) {
    if (!soft_equiv(xi[d], xi_check[d]))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("All ComptonData CSK values returned correctly.");
  } else {
    FAILMSG("Some ComptonData CSK values incorrect!");
  }

  // Finally, check the booleans
  if (!Cdata->is_legendre())
    ITFAILS;
  if (Cdata->is_lagrange())
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("ComptonData booleans set correctly.");
    std::cout << "Standard ComptonData container unit test PASSED."
              << std::endl;
  } else {
    FAILMSG("ComptonData booleans set incorrectly!");
    std::cout << "Standard ComptonData container unit test FAILED."
              << std::endl;
  }
}

void tst_lagrange_data(rtt_dsxx::ScalarUnitTest &ut) {
  std::cout << "==============================================" << std::endl;
  std::cout << "=== Testing Lagrange ComptonData container ===" << std::endl;
  std::cout << "==============================================" << std::endl;
  // make up some fake data...
  // (we tested get() and set() for the eval points in the standard test,
  // so here we only put data in the breakpoint vectors)
  std::vector<double> etemp_bp(3, 0.0);
  etemp_bp[0] = 0.3333;
  etemp_bp[1] = 0.6666;
  etemp_bp[2] = 1.9999;
  std::vector<double> gin_bp(3, 0.0);
  gin_bp[0] = 0.022;
  gin_bp[1] = 0.5555;
  gin_bp[2] = 0.988;
  std::vector<double> gout_bp(3, 0.0);
  gout_bp[0] = 0.011;
  gout_bp[1] = 0.4444;
  gout_bp[2] = 1.0;
  std::vector<double> etemp(2, 0.0);
  std::vector<double> gin(3, 0.0);
  std::vector<double> gout(2, 0.0);
  std::vector<double> xi(5, 0.0);

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      2, std::vector<std::vector<std::vector<double>>>(
             3 * 3,
             std::vector<std::vector<double>>(2, std::vector<double>(5, 0.0))));

  fake_data[0][0][0][0] = 1.345e-15;
  fake_data[0][0][1][1] = 8.888e-46;
  fake_data[0][3][0][2] = 2.687e-01;
  fake_data[0][3][1][3] = 8.888e-04;
  fake_data[1][3][0][4] = 1.000;
  fake_data[1][8][1][0] = 3.141e-29;
  fake_data[1][8][0][1] = 9.999e-99;
  fake_data[1][8][1][2] = 4.646e-02;

  // set interpolation type bools:
  bool lagrange = true;

  // try creating a compton data object:
  SP<ComptonData> Cdata;
  try {
    Cdata.reset(new ComptonData(etemp.size(), gin.size(), gout.size(),
                                xi.size(), lagrange));
    PASSMSG("ComptonData object successfully created!");
  } catch (rtt_dsxx::assertion &as) {
    FAILMSG(as.what());
    std::ostringstream message;
    message << "Aborting tests because unable to instantiate "
            << "ComptonData object";
    FAILMSG(message.str());
    return;
  }

  // next, set all of the standard and breakpoint vectors:
  Cdata->set_breakpts(etemp_bp, gin_bp, gout_bp);
  Cdata->set_evalpts(etemp, gin, gout, xi);
  Cdata->set_csk_data(fake_data);

  std::vector<double> etemp_bpcheck = Cdata->get_etemp_breakpts();
  std::vector<double> gin_bpcheck = Cdata->get_gin_breakpts();
  std::vector<double> gout_bpcheck = Cdata->get_gout_breakpts();
  std::vector<double> etemp_check = Cdata->get_etemp_pts();
  std::vector<double> gin_check = Cdata->get_gin_pts();
  std::vector<double> gout_check = Cdata->get_gout_pts();
  std::vector<double> xi_check = Cdata->get_xi_pts();
  std::vector<std::vector<std::vector<std::vector<double>>>> data_check =
      Cdata->get_csk_data();

  // check the sizes
  if (etemp_bp.size() != etemp_bpcheck.size())
    ITFAILS;
  if (gin_bp.size() != gin_bpcheck.size())
    ITFAILS;
  if (gout_bp.size() != gout_bpcheck.size())
    ITFAILS;
  if (etemp.size() != etemp_check.size())
    ITFAILS;
  if (gin.size() != gin_check.size())
    ITFAILS;
  if (gout.size() != gout_check.size())
    ITFAILS;
  if (xi.size() != xi_check.size())
    ITFAILS;
  if (Cdata->get_n_etemp_breakpts() != etemp_bpcheck.size())
    ITFAILS;
  if (Cdata->get_n_gin_breakpts() != gin_bpcheck.size())
    ITFAILS;
  if (Cdata->get_n_gout_breakpts() != gout_bpcheck.size())
    ITFAILS;
  if (Cdata->get_n_etemp_pts() != etemp_check.size())
    ITFAILS;
  if (Cdata->get_n_gin_pts() != gin_check.size())
    ITFAILS;
  if (Cdata->get_n_gout_pts() != gout_check.size())
    ITFAILS;
  if (Cdata->get_n_xi_pts() != xi_check.size())
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
    for (size_t b = 0; b < 3 * gin.size(); b++) {
      for (size_t c = 0; c < gout.size(); c++) {
        for (size_t d = 0; d < xi.size(); d++) {
          if (!soft_equiv(data_check[a][b][c][d], fake_data[a][b][c][d]))
            ITFAILS;
        }
      }
    }
  }

  // check etemp, gin, and gout breakpoints:
  for (size_t b = 0; b < etemp_bp.size(); b++) {
    if (!soft_equiv(etemp_bp[b], etemp_bpcheck[b]))
      ITFAILS;
  }
  for (size_t b = 0; b < gin_bp.size(); b++) {
    if (!soft_equiv(gin_bp[b], gin_bpcheck[b]))
      ITFAILS;
  }
  for (size_t c = 0; c < gout_bp.size(); c++) {
    if (!soft_equiv(gout_bp[c], gout_bpcheck[c]))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("All ComptonData CSK values returned correctly.");
  } else {
    FAILMSG("Some ComptonData CSK values incorrect!");
  }

  // Finally, check the booleans
  if (Cdata->is_legendre())
    ITFAILS;
  if (!Cdata->is_lagrange())
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("ComptonData booleans set correctly.");
    std::cout << "Lagrange ComptonData container unit test PASSED."
              << std::endl;
  } else {
    FAILMSG("ComptonData booleans set incorrectly!");
    std::cout << "Lagrange ComptonData container unit test FAILED."
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
  std::vector<double> gin_bp(3, 0.0);
  gin_bp[0] = 0.022;
  gin_bp[1] = 0.5555;
  gin_bp[2] = 0.988;
  std::vector<double> gout_bp(3, 0.0);
  gout_bp[0] = 0.011;
  gout_bp[1] = 0.4444;
  gout_bp[2] = 1.0;
  std::vector<double> etemp(2, 0.0);
  etemp[0] = 0.4444;
  etemp[1] = 0.7777;
  std::vector<double> gin(4, 0.0);
  gin[0] = 0.25;
  gin[1] = 0.35;
  gin[2] = 0.65;
  gin[3] = 0.75;
  std::vector<double> gout(4, 0.0);
  gout[0] = 0.2;
  gout[1] = 0.3;
  gout[2] = 0.6;
  gout[3] = 0.7;
  std::vector<double> xi(1, 0.0);
  xi[0] = 0.0;

  // Fill in certain points of the fake data...
  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
      2, std::vector<std::vector<std::vector<double>>>(
             4 * 3,
             std::vector<std::vector<double>>(4, std::vector<double>(1, 0.0))));

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

  // set interpolation type bools:
  bool lagrange = true;
  bool legendre = false;

  // try creating a compton data object:
  SP<ComptonData> Cdata;
  try {
    Cdata.reset(new ComptonData(etemp.size(), gin.size(), gout.size(),
                                xi.size(), lagrange, legendre));
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
  Cdata->set_breakpts(etemp_bp, gin_bp, gout_bp);
  Cdata->set_evalpts(etemp, gin, gout, xi);
  Cdata->set_csk_data(fake_data);

  // get the slice of data corresponding to electron temperature index 0:
  std::vector<std::vector<std::vector<std::vector<double>>>> slice0 =
      Cdata->get_csk_data(0);

  // get the slice of data corresponding to electron temperature index 1:
  std::vector<std::vector<std::vector<std::vector<double>>>> slice1 =
      Cdata->get_csk_data(1);

  // check the returned data against the original values
  for (size_t a = 0; a < slice1.size(); a++) {
    for (size_t b = 0; b < slice1[a].size(); b++) {
      for (size_t c = 0; c < slice1[a][b].size(); c++) {
        for (size_t d = 0; d < slice1[a][b][c].size(); d++) {
          if (!soft_equiv(slice0[a][b][c][d], fake_data[0][b][c][d]))
            ITFAILS;
          if (!soft_equiv(slice1[a][b][c][d], fake_data[1][b][c][d]))
            ITFAILS;
        }
      }
    }
  }
  if (ut.numFails == 0) {
    PASSMSG("All ComptonData large data slices returned correctly.");
  } else {
    FAILMSG("Some ComptonData large data slice values incorrect!");
  }

  // get the slice of data corresponding to electron temperature index 0,
  // interpolation region 0 (bottom), x-region 1, y-region 1:
  std::vector<std::vector<std::vector<std::vector<double>>>> slice0011 =
      Cdata->get_csk_data(0, 0, 1, 1);

  // get the slice of data corresponding to electron temperature index 1,
  // interpolation region 1 (middle), x-region 0, y-region 0:
  std::vector<std::vector<std::vector<std::vector<double>>>> slice1100 =
      Cdata->get_csk_data(1, 1, 0, 0);

  // check slice size-- should be ppr_etemp x ppr_gin x ppr_gout x n_xi
  // (which is 1 x 2 x 2 x 1)
  size_t ppr_gout = gout.size() / (gout_bp.size() - 1);
  size_t ppr_gin = gout.size() / (gout_bp.size() - 1);

  // check the returned data against the original values
  for (size_t a = 0; a < slice0011.size(); a++) {
    for (size_t b = 0; b < slice0011[a].size(); b++) {
      for (size_t c = 0; c < slice0011[a][b].size(); c++) {
        for (size_t d = 0; d < slice0011[a][b][c].size(); d++) {
          if (!soft_equiv(slice0011[a][b][c][d],
                          fake_data[0][ppr_gin + b][ppr_gout + c][d]))
            ITFAILS;
        }
      }
    }
  }

  // check the returned data against the original values
  for (size_t a = 0; a < slice1100.size(); a++) {
    for (size_t b = 0; b < slice1100[a].size(); b++) {
      for (size_t c = 0; c < slice1100[a][b].size(); c++) {
        for (size_t d = 0; d < slice1100[a][b][c].size(); d++) {
          if (!soft_equiv(slice1100[a][b][c][d],
                          fake_data[1][2 * ppr_gin + b][c][d]))
            ITFAILS;
        }
      }
    }
  }
  if (ut.numFails == 0) {
    PASSMSG("All ComptonData single-region data slices returned correctly.");
  } else {
    FAILMSG("Some ComptonData single-region data slice values incorrect!");
  }
}

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_accessors(ut);
    tst_std_data(ut);
    tst_lagrange_data(ut);
    tst_data_slices(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
