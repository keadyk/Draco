//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/tstComptonFile.cc
 * \author Kendra P. Keady
 * \date   Mon Oct 17 2016
 * \brief  Compton test.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
#include "ComptonFile.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/SP.hh"
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

using rtt_compton::ComptonFile;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

// prototype functions to check data:
void check_points(rtt_dsxx::ScalarUnitTest&, const std::vector<double>&,
                  const std::vector<double>&, const std::vector<double>&, 
                  const std::vector<double>&);

void check_data(rtt_dsxx::ScalarUnitTest&, 
             const std::vector<std::vector<std::vector<std::vector<double>>>>&);

void tst_asciiread(rtt_dsxx::ScalarUnitTest& ut)
{
  std::cout << "======= Testing ascii CSK library read =======" << std::endl;
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
  std::vector<std::vector<std::vector<std::vector<double>>>>tst_data = 
    spCompton->read_csk_data();

  std::vector<double> gin_pts = spCompton->get_gin_pts();
  std::vector<double> gout_pts = spCompton->get_gout_pts();
  std::vector<double> xi_pts = spCompton->get_xi_pts();
  std::vector<double> etemp_pts = spCompton->get_etemp_pts();

  // check the points to be sure they match the file:
  check_points(ut, gin_pts, gout_pts, xi_pts, etemp_pts);
  check_data(ut, tst_data);

}

void tst_binread(rtt_dsxx::ScalarUnitTest& ut)
{
  std::cout << "======= Testing binary CSK library read =======" << std::endl;
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
  std::vector<std::vector<std::vector<std::vector<double>>>>tst_data = 
    spCompton->read_csk_data();

  std::vector<double> gin_pts = spCompton->get_gin_pts();
  std::vector<double> gout_pts = spCompton->get_gout_pts();
  std::vector<double> xi_pts = spCompton->get_xi_pts();
  std::vector<double> etemp_pts = spCompton->get_etemp_pts();

  // check the points to be sure they match the file:
  check_points(ut, gin_pts, gout_pts, xi_pts, etemp_pts);
  check_data(ut, tst_data);

  
}

void check_points(rtt_dsxx::ScalarUnitTest& ut, 
                    const std::vector<double>&gin_pts,
                    const std::vector<double>&gout_pts, 
                    const std::vector<double>&xi_pts,
                    const std::vector<double>&etemp_pts)
{
  double gin_ref, etemp_ref;
  // check the uniformly-discretized variables:
  for(size_t a=0; a<20; a++)
  {
    gin_ref = 0.05 + a*0.1;
    if(!soft_equiv(gin_pts[a], gin_ref)) ITFAILS;;
  }
  // finely-discretized region
  for(size_t b=0; b<5; b++)
  {
    etemp_ref = 0.025 + b*0.05;
    if(!soft_equiv(etemp_pts[b], etemp_ref)) ITFAILS;;
  }  
  // coarsely-discretized region
  for(size_t c=0; c<5; c++)
  {
    etemp_ref = 0.325 + c*0.15;
    if(!soft_equiv(etemp_pts[5+c], etemp_ref)) ITFAILS;;
  }  

  // check the arbitrarily-discretized data
  if(!soft_equiv(gout_pts[0], 0.01)) ITFAILS;
  if(!soft_equiv(gout_pts[1], 0.035)) ITFAILS;
  if(!soft_equiv(gout_pts[2], 0.28)) ITFAILS;
  if(!soft_equiv(gout_pts[3], 0.5)) ITFAILS;
  if(!soft_equiv(gout_pts[4], 0.99)) ITFAILS;

  if(!soft_equiv(xi_pts[0], -0.95)) ITFAILS;
  if(!soft_equiv(xi_pts[1], 0.0)) ITFAILS;
  if(!soft_equiv(xi_pts[2], 0.25)) ITFAILS;
  if(!soft_equiv(xi_pts[3], 0.88)) ITFAILS;
  if(!soft_equiv(xi_pts[4], 0.95)) ITFAILS;

  if (ut.numFails==0)
  {
    PASSMSG("Successfully read evaluation points from test library!");
  }
  else
  {
    FAILMSG("Did not read correct evaluation points from test library!");
  } 
}

void check_data(rtt_dsxx::ScalarUnitTest& ut,
         const std::vector<std::vector<std::vector<std::vector<double>>>>& data)
{
  // check selected data points from each regime:
  if(data.size() != 10) ITFAILS;
  if(data[0].size() != 20) ITFAILS;
  if(data[0][0].size() != 5) ITFAILS;
  if(data[0][0][0].size() != 5) ITFAILS;
  // SCALE small values and compare to one:
  if(!soft_equiv(data[0][0][0][0], 5.55794e-05, 1e-5)) ITFAILS;  
  if(!soft_equiv(data[0][0][0][4]/2.02745e-81, 1.0, 1e-5)) ITFAILS; 
  if(!soft_equiv(data[0][0][4][0]/3.07426e-31, 1.0, 1e-5)) ITFAILS;
  if(!soft_equiv(data[0][0][4][4]/2.9701e-222, 1.0, 1e-5)) ITFAILS;
  if(!soft_equiv(data[0][19][0][0]/1.07261e-92, 1.0, 1e-5)) ITFAILS;
  if(!soft_equiv(data[0][19][0][4], 0.0, 1e-5)) ITFAILS;
  if(!soft_equiv(data[0][19][4][0], 2.0226e-05, 1e-5)) ITFAILS;
  if(!soft_equiv(data[0][19][4][4]/3.45633e-17, 1.0, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][0][0][0], 0.755117, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][0][0][4], 0.00901667, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][0][4][0], 0.862336, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][0][4][4], 6.41622e-06, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][19][0][0], 3.62824e-06, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][19][0][4]/2.31686e-23, 1.0, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][19][4][0], 0.124555, 1e-5)) ITFAILS;
  if(!soft_equiv(data[9][19][4][4], 0.0684458, 1e-5)) ITFAILS;

  if (ut.numFails==0)
  {
    PASSMSG("Successfully read CSK data points from test library!");
  }
  else
  {
    FAILMSG("Did not read correct CSK data points from test library!");
  } 
}

int main(int argc, char *argv[])
{
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    tst_asciiread(ut);
    tst_binread(ut);
  }
  UT_EPILOG(ut);
  return 0;
}
