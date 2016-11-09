//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/tstComptonInterp.cc
 * \author Kendra P. Keady
 * \date   Mon Oct 17 2016
 * \brief  Compton interpolation test.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */

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
