//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/tstComptonInterp.cc
 * \author Kendra P. Keady
 * \date   Mon Oct 17 2016
 * \brief  Compton interpolation test.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */

#include "ComptonData.hh"
#include "ComptonInterp.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/SP.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cstring>
#include <iostream>
#include <vector>

using rtt_compton::ComptonInterp;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

// This function tests the one-dimensional (electron temperature)
// Lagrange interpolation capability in Draco. We use the fact that Lagrange
// interpolation with N points interpolates an (N-1)-degree polynomial exactly
// to check the implementation.
void tst_1D_interp(rtt_dsxx::ScalarUnitTest &ut) {
  // TEST ONE: Interpolation of linear CSK data.
  {
    std::cout << "==============================================" << std::endl;
    std::cout << "= Testing electron temperature interpolation =" << std::endl;
    std::cout << "==============================================" << std::endl;
    // make a phony comptondata object and fill in part of its data
    size_t ngin = 2;
    size_t ngout = 2;
    size_t netemp = 2;
    size_t nxi = 1;

    // set only one region per variable:
    size_t nginbp = 2;
    size_t ngoutbp = 2;
    size_t netempbp = 2;

    bool lag = true;
    bool leg = false;
    rtt_dsxx::SP<rtt_compton::ComptonData> Cdata(
        new rtt_compton::ComptonData(netemp, ngin, ngout, nxi, lag, leg));
    std::vector<double> ginpts(ngin, 0.0);
    std::vector<double> goutpts(ngout, 0.0);
    std::vector<double> etemppts(netemp, 0.0);
    std::vector<double> xipts(nxi, 0.0);

    std::vector<double> ginbpts(nginbp, 0.0);
    std::vector<double> goutbpts(ngoutbp, 0.0);
    std::vector<double> etempbpts(netempbp, 0.0);

    // fill in the data:
    ginpts[0] = 0.25;
    ginpts[1] = 0.75;

    goutpts[0] = 0.25;
    goutpts[1] = 0.75;

    etemppts[0] = 0.1;
    etemppts[1] = 0.9;

    xipts[0] = 0.0;

    ginbpts[0] = 0.0;
    ginbpts[1] = 1.0;

    goutbpts[0] = 0.0;
    goutbpts[1] = 1.0;

    etempbpts[0] = 0.0;
    etempbpts[1] = 1.0;

    // set all of this data:
    Cdata->set_evalpts(etemppts, ginpts, goutpts, xipts);
    Cdata->set_breakpts(etempbpts, ginbpts, goutbpts);

    std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
        netemp, std::vector<std::vector<std::vector<double>>>(
                    3 * ngin, std::vector<std::vector<double>>(
                                  ngout, std::vector<double>(nxi, 0.0))));

    // now, form the actual csk data:
    for (size_t a = 0; a < netemp; a++) {
      for (size_t b = 0; b < 3 * ngin; b++) {
        for (size_t c = 0; c < ngout; c++) {
          for (size_t d = 0; d < nxi; d++) {
            fake_data[a][b][c][d] =
                static_cast<double>(a * (3 * 4) + b * 4 + c) + 0.5;
          }
        }
      }
    }
    Cdata->set_csk_data(fake_data);

    // try to make a compton interpolation object
    ComptonInterp Cinterp(Cdata);

    // pick an electron temperature w/in the breakpoint bounds:
    double test_etemp1 = 0.55;

    // interpolate ALL data in gin, gout, and xi for the given etemp:
    std::vector<std::vector<std::vector<double>>> interp_etemp1 =
        Cinterp.interpolate_etemp(test_etemp1);

    // Check the results:
    for (size_t g = 0; g < interp_etemp1.size(); g++) {
      for (size_t h = 0; h < interp_etemp1[g].size(); h++) {
        double lin_value = fake_data[0][g][h][0] +
                           (0.55 - etemppts[0]) / (etemppts[1] - etemppts[0]) *
                               (fake_data[1][g][h][0] - fake_data[0][g][h][0]);
        if (!soft_equiv(interp_etemp1[g][h][0], lin_value))
          ITFAILS;
      }
    }

    if (ut.numFails == 0) {
      PASSMSG("ComptonInterp linear function interpolation okay.");
    } else {
      FAILMSG("ComptonInterp linear function interpolation failed!");
    }
  }

  // TEST TWO: Interpolation of fourth-degree polynomial CSK data--
  // requires five electron temperature points to be exact.
  {
    // make a phony comptondata object and fill in part of its data
    size_t ngin = 1;
    size_t ngout = 1;
    size_t netemp = 5;
    size_t nxi = 1;

    // set only one region per variable:
    size_t nginbp = 2;
    size_t ngoutbp = 2;
    size_t netempbp = 2;

    bool lag = true;
    bool leg = false;
    rtt_dsxx::SP<rtt_compton::ComptonData> Cdata(
        new rtt_compton::ComptonData(netemp, ngin, ngout, nxi, lag, leg));
    std::vector<double> ginpts(ngin, 0.0);
    std::vector<double> goutpts(ngout, 0.0);
    std::vector<double> etemppts(netemp, 0.0);
    std::vector<double> xipts(nxi, 0.0);

    std::vector<double> ginbpts(nginbp, 0.0);
    std::vector<double> goutbpts(ngoutbp, 0.0);
    std::vector<double> etempbpts(netempbp, 0.0);

    // fill in the data:
    ginpts[0] = 0.5;
    goutpts[0] = 0.5;

    etemppts[0] = 0.1;
    etemppts[1] = 0.2;
    etemppts[2] = 0.6;
    etemppts[3] = 0.75;
    etemppts[4] = 0.9;

    xipts[0] = 0.0;

    ginbpts[0] = 0.0;
    ginbpts[1] = 1.0;

    goutbpts[0] = 0.0;
    goutbpts[1] = 1.0;

    etempbpts[0] = 0.0;
    etempbpts[1] = 1.0;

    // set all of this data:
    Cdata->set_evalpts(etemppts, ginpts, goutpts, xipts);
    Cdata->set_breakpts(etempbpts, ginbpts, goutbpts);

    std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
        netemp, std::vector<std::vector<std::vector<double>>>(
                    3 * ngin, std::vector<std::vector<double>>(
                                  ngout, std::vector<double>(nxi, 0.0))));

    // now, form the actual csk data, choosing the points to fit a fourth-order
    // polynomial
    for (size_t a = 0; a < netemp; a++) {
      double epow4 = std::pow(etemppts[a], 4);
      double epow3 = std::pow(etemppts[a], 3);
      double epow2 = etemppts[a] * etemppts[a];
      fake_data[a][0][0][0] = 2.0 * epow4 + 0.3 * epow3 + 0.484;
      fake_data[a][1][0][0] =
          1.5 * epow4 + 1.4 * epow2 - 0.2 * etemppts[a] + 0.333;
      fake_data[a][2][0][0] = 0.8 * epow4;
    }
    Cdata->set_csk_data(fake_data);

    // try to make a compton interpolation object
    ComptonInterp Cinterp(Cdata);

    // pick two electron temperature w/in the breakpoint bounds:
    double test_etemp1 = 0.35;
    double test_etemp2 = 0.85;

    // interpolate ALL data in gin, gout, and xi for the given etemp:
    std::vector<std::vector<std::vector<double>>> interp_etemp1 =
        Cinterp.interpolate_etemp(test_etemp1);
    std::vector<std::vector<std::vector<double>>> interp_etemp2 =
        Cinterp.interpolate_etemp(test_etemp2);

    // Check the results for etemp=0.35:
    double poly_value11 =
        2.0 * std::pow(test_etemp1, 4) + 0.3 * std::pow(test_etemp1, 3) + 0.484;
    double poly_value12 = 1.5 * std::pow(test_etemp1, 4) +
                          1.4 * std::pow(test_etemp1, 2) - 0.2 * test_etemp1 +
                          0.333;
    double poly_value13 = 0.8 * std::pow(test_etemp1, 4);

    if (!soft_equiv(interp_etemp1[0][0][0], poly_value11))
      ITFAILS;
    if (!soft_equiv(interp_etemp1[1][0][0], poly_value12))
      ITFAILS;
    if (!soft_equiv(interp_etemp1[2][0][0], poly_value13))
      ITFAILS;

    // Check the results for etemp=0.85:
    double poly_value21 =
        2.0 * std::pow(test_etemp2, 4) + 0.3 * std::pow(test_etemp2, 3) + 0.484;
    double poly_value22 = 1.5 * std::pow(test_etemp2, 4) +
                          1.4 * std::pow(test_etemp2, 2) - 0.2 * test_etemp2 +
                          0.333;
    double poly_value23 = 0.8 * std::pow(test_etemp2, 4);

    if (!soft_equiv(interp_etemp2[0][0][0], poly_value21))
      ITFAILS;
    if (!soft_equiv(interp_etemp2[1][0][0], poly_value22))
      ITFAILS;
    if (!soft_equiv(interp_etemp2[2][0][0], poly_value23))
      ITFAILS;
  }

  if (ut.numFails == 0) {
    PASSMSG("ComptonInterp fourth-order polynomial interpolation okay.");
    std::cout << "Electron temperature interpolation PASSED." << std::endl;
  } else {
    FAILMSG("ComptonInterp fourth-order polynomial interpolation failed!");
    std::cout << "Electron temperature interpolation FAILED." << std::endl;
  }
}

// This function tests the two-dimensional (gamma_in and gamma_out)
// Lagrange interpolation capability in Draco. We use the fact that Lagrange
// interpolation with N points interpolates an (N-1)-degree polynomial exactly
// to check the implementation.
void tst_2D_interp(rtt_dsxx::ScalarUnitTest &ut) {
  // TEST ONE: Interpolation of linear CSK data.
  {
    std::cout << "==============================================" << std::endl;
    std::cout << "=== Testing frequency in/out interpolation ===" << std::endl;
    std::cout << "==============================================" << std::endl;
    // make a phony comptondata object and fill in part of its data
    size_t ngin = 2;
    size_t ngout = 2;
    size_t netemp = 1;
    size_t nxi = 1;

    // set only one region per variable:
    size_t nginbp = 2;
    size_t ngoutbp = 2;
    size_t netempbp = 2;

    bool lag = true;
    bool leg = false;
    rtt_dsxx::SP<rtt_compton::ComptonData> Cdata(
        new rtt_compton::ComptonData(netemp, ngin, ngout, nxi, lag, leg));
    std::vector<double> ginpts(ngin, 0.0);
    std::vector<double> goutpts(ngout, 0.0);
    std::vector<double> etemppts(netemp, 0.0);
    std::vector<double> xipts(nxi, 0.0);

    std::vector<double> ginbpts(nginbp, 0.0);
    std::vector<double> goutbpts(ngoutbp, 0.0);
    std::vector<double> etempbpts(netempbp, 0.0);

    // fill in the data:
    ginpts[0] = 0.20;
    ginpts[1] = 0.80;
    goutpts[0] = 0.20;
    goutpts[1] = 0.80;

    etemppts[0] = 0.5;

    xipts[0] = 0.0;

    ginbpts[0] = 0.0;
    ginbpts[1] = 1.0;

    goutbpts[0] = 0.0;
    goutbpts[1] = 1.0;

    etempbpts[0] = 0.0;
    etempbpts[1] = 1.0;

    // set all of this data:
    Cdata->set_evalpts(etemppts, ginpts, goutpts, xipts);
    Cdata->set_breakpts(etempbpts, ginbpts, goutbpts);

    std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(
        netemp, std::vector<std::vector<std::vector<double>>>(
                    3 * ngin, std::vector<std::vector<double>>(
                                  ngout, std::vector<double>(nxi, 0.0))));

    // now, form the actual csk data:
    // (take advantage of the fact that the vector is size 1 in xi and etemp
    // to simplify the loops)
    for (size_t b = 0; b < 3 * ngin; b++) {
      fake_data[0][b][0][0] = 1.25 * (b + 0.5);
      fake_data[0][b][1][0] = 1.25 * (b + 0.5) + 0.5;
    }

    Cdata->set_csk_data(fake_data);

    // grab the etemp-0 data:
    std::vector<std::vector<std::vector<double>>> input_data = fake_data[0];

    // try to make a compton interpolation object
    ComptonInterp Cinterp(Cdata);

    // pick gamma in/out values to test all three boundary-layer regions:
    double test_gin1 = 0.33;
    double test_gin2 = 0.6;
    double test_gin3 = 0.5;
    double test_gin4 = 0.25;

    double test_gout1 = 0.1;
    double test_gout2 = 0.4;
    double test_gout3 = 0.52;
    double test_gout4 = 0.77;

    // interpolate ALL data in xi for the given gin/gout pair
    // (just 1 point in xi for the case considered here):
    std::vector<double> interp_ginout1 =
        Cinterp.interpolate_gin_gout(test_gin1, test_gout1, input_data);
    std::vector<double> interp_ginout2 =
        Cinterp.interpolate_gin_gout(test_gin2, test_gout1, input_data);
    std::vector<double> interp_ginout3 =
        Cinterp.interpolate_gin_gout(test_gin2, test_gout2, input_data);
    std::vector<double> interp_ginout4 =
        Cinterp.interpolate_gin_gout(test_gin3, test_gout2, input_data);
    std::vector<double> interp_ginout5 =
        Cinterp.interpolate_gin_gout(test_gin3, test_gout3, input_data);
    std::vector<double> interp_ginout6 =
        Cinterp.interpolate_gin_gout(test_gin4, test_gout4, input_data);

    // while we're at it, check to see that the original x/y points correctly
    // evaluate to the CSK data... (or get close, at least)
    std::vector<double> interp_ginout7 =
        Cinterp.interpolate_gin_gout(0.20000001, 1.0000001, input_data);
    std::vector<double> interp_ginout8 =
        Cinterp.interpolate_gin_gout(0.2000001, 0.4000001, input_data);

    // compare the interp results to the reference cases for the low,
    // middle, and high regions:
    double x_coeff = 2.083333333;
    double y_coeff = 0.833333333;
    double c_lo = 0.041666667; //(green. c lo green.)
    double c_mid = 2.541666667;
    double c_hi = 5.041666667;

    // calculate the reference values (remembering to conver gin/out into x/y):
    double val1 =
        x_coeff * test_gin1 +
        y_coeff *
            (1.0 - (2.0 * test_gin1 * test_gout1) / (test_gin1 - test_gout1)) +
        c_lo;
    double val2 =
        x_coeff * test_gin2 +
        y_coeff *
            (1.0 - (2.0 * test_gin2 * test_gout1) / (test_gin2 - test_gout1)) +
        c_lo;
    double val3 =
        x_coeff * test_gin2 +
        y_coeff *
            (1.0 - (test_gin2 - test_gout2) / (2.0 * test_gin2 * test_gout2)) +
        c_mid;
    double val4 =
        x_coeff * test_gin3 +
        y_coeff *
            (1.0 - (test_gin3 - test_gout2) / (2.0 * test_gin3 * test_gout2)) +
        c_mid;
    double val5 =
        x_coeff * test_gin3 + y_coeff * (test_gout3 - test_gin3) + c_hi;
    double val6 =
        x_coeff * test_gin4 + y_coeff * (test_gout4 - test_gin4) + c_hi;

    if (!soft_equiv(val1, interp_ginout1[0], 1e-7)) {
      std::cout << val1 << " != " << interp_ginout1[0] << std::endl;
      FAILMSG("Interp failed in bottom boundary-layer region!");
    }
    if (!soft_equiv(val2, interp_ginout2[0], 1e-7)) {
      std::cout << val2 << " != " << interp_ginout2[0] << std::endl;
      FAILMSG("Interp failed in bottom boundary-layer region!");
    }
    if (!soft_equiv(val3, interp_ginout3[0], 1e-7)) {
      std::cout << val3 << " != " << interp_ginout3[0] << std::endl;
      FAILMSG("Interp failed in middle boundary-layer region!");
    }
    if (!soft_equiv(val4, interp_ginout4[0], 1e-7)) {
      std::cout << val4 << " != " << interp_ginout4[0] << std::endl;
      FAILMSG("Interp failed in middle boundary-layer region!");
    }
    if (!soft_equiv(val5, interp_ginout5[0], 1e-7)) {
      std::cout << val5 << " != " << interp_ginout5[0] << std::endl;
      FAILMSG("Interp failed in top boundary-layer region!");
    }
    if (!soft_equiv(val6, interp_ginout6[0], 1e-7)) {
      std::cout << val6 << " != " << interp_ginout6[0] << std::endl;
      FAILMSG("Interp failed in top boundary-layer region!");
    }
    if (!soft_equiv(6.125, interp_ginout7[0], 1e-4)) {
      std::cout << 6.125 << " != " << interp_ginout7[0] << std::endl;
      FAILMSG("Interp failed to reproduce CSK data point!");
    }
    if (!soft_equiv(5.625, interp_ginout8[0], 1e-4)) {
      std::cout << 5.625 << " != " << interp_ginout8[0] << std::endl;
      FAILMSG("Interp failed to reproduce CSK data point!");
    }

    if (ut.numFails == 0) {
      PASSMSG("ComptonInterp 2-D linear function interpolation okay.");
    } else {
      FAILMSG("ComptonInterp 2-D linear function interpolation failed!");
    }
  }
}

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_1D_interp(ut);
    tst_2D_interp(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
