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
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <string>
#include <iostream>

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
 * \arg interpolation data?
 *
 */

namespace rtt_compton {

// Enumerated type to describe the interpolation regions
// (used only within this class)
enum class Region{ BOTTOM=1, MID=2, TOP=3, NONE=4 };

class ComptonInterp {
  private:
    //! smart pointer to constant compton data:
    rtt_dsxx::SP<const ComptonData> Cdata;

    //! function to determine what gamma/gammaout region we're in
    //Region find_global_region(const double, const double);

    //! function to determine what electron temp interpolation region we're in
    //double find_local_region(const double);

    //! function to determine what local x/y interpolation region we're in
    //std::pair<double, double> find_local_region(const double, const double);
  public:
    // Constructor 
    // TODO: I think all of the vectors of data should be placed in a compton
    // data class, which can then simply be passed to this class via a SP.
    ComptonInterp(rtt_dsxx::SP<const ComptonData>);

    // Destructor
    ~ComptonInterp();

    //! Interpolate in electron temperature for a CSK point (or integral)
    /*double interpolate_etemp(double, double, double);

    //! Interpolate in electron temperature for all CSK points in gin/gout
    std::vector<double> interpolate_etemp(std::vector<double>, 
                                          std::vector<double>,
                                          double);*/
};

}
#endif
