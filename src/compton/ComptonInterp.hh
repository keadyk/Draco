//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/ComptonInterp.hh
 * \author Kendra Keady
 * \date   Tues Nov 8 2016 (Election Day!)
 * \brief  Header file for ComptonInterp class
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
#ifndef __compton_ComptonInterp_hh__
#define __compton_ComptonInterp_hh__

#include "ComptonData.hh"
#include "ds++/SP.hh"
#include <cmath>
#include <iostream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

//===========================================================================//
/*!
 * \class ComptonInterp
 *
 * \brief Class to interpolate between Compton Scattering Kernel (CSK) values
 *
 * The ComptonInterp class provides routines to interpolate between adjacent
 * CSK values in temperature, incident frequency, and outgoing frequency.
 * 
 *
 * \arg SP<ComptonData> - smart pointer to a Compton Data object
 *
 */

namespace rtt_compton {

// Enumerated type to describe the interpolation regions
// (used only within this class)
enum class Region { BOTTOM = 0, MID = 1, TOP = 2, NONE = 3 };

class ComptonInterp {
private:
  //! smart pointer to constant compton data:
  rtt_dsxx::SP<const ComptonData> Cdata;

  //! Product values for Lagrange interpolation (in x, y. and etemp)
  std::vector<double> prod_x, prod_y, prod_etemp;

  //! scaled and shifted x, y, and etemp values:
  std::vector<double> xs, ys, etemps;

  //! Number of breakpoints per variable:
  size_t nx_break, ny_break, netemp_break;
  //! Number of local interpolation points per breakpoint region:
  size_t nx_local, ny_local, netemp_local;

  size_t binary_search(const double, const std::vector<double> &);

  Region find_global_region(const double, const double);

  //! function to determine what electron temp interpolation region we're in
  size_t find_etemp_region(const double);

  //! function to determine what local x/y interpolation region we're in
  void set_xy_and_region(const double, const double, const Region,
                         std::pair<size_t, size_t> &,
                         std::pair<double, double> &);

public:
  // Constructor
  ComptonInterp(rtt_dsxx::SP<const ComptonData>);

  // Destructor
  ~ComptonInterp(){};

  //! Interpolate ALL gin/gout/xi data in electron temperature
  std::vector<std::vector<std::vector<double>>> interpolate_etemp(const double);

  //! Interpolate ALL xi data for a gin/gout
  std::vector<double>
  interpolate_gin_gout(const double, const double,
                       const std::vector<std::vector<std::vector<double>>> &);

  //! Use Lagrange interpolation on a single etemp point
  double interpolate_etemp(const std::vector<double> &,
                           const std::vector<double> &, const double);

  //! Use Lagrange interpolation on a single frequency in/out pair
  double interpolate_gin_gout(const std::vector<std::vector<double>> &,
                              const std::vector<double> &,
                              const std::vector<double> &, const double,
                              const double);
};
}
#endif
