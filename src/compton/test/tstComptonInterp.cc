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
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/SP.hh"
#include <cstring>
#include <iostream>
#include <vector>

using rtt_compton::ComptonInterp;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

void tst_interp(rtt_dsxx::ScalarUnitTest& ut)
{
  {
  // make a phony comptondata object and fill in part of its data
  size_t ngin = 2;
  size_t ngout = 2;
  size_t netemp = 5;
  size_t nxi = 1;

  // set only one region per variable:
  size_t nginbp = 2;
  size_t ngoutbp = 2;
  size_t netempbp = 2;

  bool lag = true;
  bool leg = false;
  rtt_dsxx::SP<rtt_compton::ComptonData> Cdata( new rtt_compton::ComptonData(netemp, ngin, ngout, nxi, lag, leg));
  std::vector<double>ginpts(ngin, 0.0);
  std::vector<double>goutpts(ngout, 0.0);
  std::vector<double>etemppts(netemp, 0.0);
  std::vector<double>xipts(nxi, 0.0);

  std::vector<double>ginbpts(nginbp, 0.0);
  std::vector<double>goutbpts(ngoutbp, 0.0);
  std::vector<double>etempbpts(netempbp, 0.0);
  
  // fill in the data:
  ginpts[0] = 0.25;
  ginpts[1] = 0.75;

  goutpts[0] = 0.25;
  goutpts[1] = 0.75;

  etemppts[0] = 0.1;
  etemppts[1] = 0.3;
  etemppts[2] = 0.5;
  etemppts[3] = 0.7;
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

  std::vector<std::vector<std::vector<std::vector<double>>>> fake_data(netemp,
    std::vector<std::vector<std::vector<double>>>(ngin, 
    std::vector<std::vector<double>>(ngout, std::vector<double>(nxi, 0.0))));

  // now, form the actual csk data:
  for(size_t a = 0; a < netemp; a++)
  {
    for(size_t b = 0; b < ngin; b++)
    {
      for(size_t c = 0; c < ngout; c++)
      {
        for(size_t d = 0; d < nxi; d++)
        {
          fake_data[a][b][c][d] = static_cast<double>(a*(12*4)+b*4+c) + 0.5;
        }
      }
    }
  }
  Cdata->set_csk_data(fake_data);

  // try to make a compton interpolation object
  ComptonInterp Cinterp(Cdata);
  
  // pick an electron temperature w/in the breakpoint bounds:
  double test_etemp = 0.55;

  // interpolate ALL data in gin, gout, and xi for the given etemp:
  std::vector<std::vector<std::vector<double>>> interp_etemp1 =
    Cinterp.interpolate_etemp(test_etemp);

  }

  ITFAILS;
}

int main(int argc, char *argv[])
{
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_interp(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
