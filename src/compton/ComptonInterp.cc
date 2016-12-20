//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/ComptonInterp.cc
 * \author Kendra Keady
 * \date   Tues Nov 8 2016 (Election Day!)
 * \brief  Class file for ComptonInterp 
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */

#include "ComptonInterp.hh"
#include "ComptonData.hh"
#include <sstream>

namespace rtt_compton {

ComptonInterp::ComptonInterp(rtt_dsxx::SP<const ComptonData> Cdata_)
    : Cdata(Cdata_) {
  etemps = Cdata->get_etemp_pts();

  netemp_local = Cdata->get_n_etemp_pts() / (Cdata->get_n_etemp_breakpts() - 1);

  netemp_break = Cdata->get_n_etemp_breakpts();

  prod_etemp.assign(etemps.size(), 0.0);

  // Pre-compute prod_etemps_{i,j} = product_{l \neq j} 1/(etemps[i,j] - etemps[i,l])
  for (size_t i = 0; i < netemp_break - 1; ++i) {
    for (size_t j = 0; j < netemp_local; ++j) {
      size_t k = j + netemp_local * i;
      prod_etemp[k] = 1.;
      for (size_t l = 0; l < netemp_local; ++l) {
        if (l != j) {
          size_t k1 = l + netemp_local * i;
          prod_etemp[k] = prod_etemp[k] / (etemps[k] - etemps[k1]);
        }
      }
    }
  }
}

size_t ComptonInterp::find_etemp_region(const double etemp) {
  std::vector<double> breakpt_data = Cdata->get_etemp_breakpts();
  size_t iregion = binary_search(etemp, breakpt_data);

  return iregion;
}

// binary search for grid index
// TODO: could be templated if we need more than just the double version
size_t ComptonInterp::binary_search(const double value,
                                    const std::vector<double> &bp_data) {
  // check that the value is actually between the bounds:
  Ensure(value >= bp_data[0] && value <= bp_data.back());

  // if there are only two options, the value had better be in region 0...
  if (bp_data.size() == 2) {
    return 0;
  }

  // get low and high indices for the breakpoint vector
  size_t low = 0;
  size_t high = bp_data.size() - 1;
  // return value (actual region index)
  size_t index;

  // binary search
  while ((high - low) > 1) {
    index = (high - low) / 2;
    if (value < bp_data[index]) {
      high = index;
    } else {
      low = index;
    }
  }

  // since we're looking for a bin index, we return the index of the lower edge.
  index = low;
  return index;
}

// interpolate ALL gin/gout/xi data for an electron temperature
std::vector<std::vector<std::vector<double>>>
ComptonInterp::interpolate_etemp(const double etemp) {
  // figure out what interpolation region we're in
  size_t i = find_etemp_region(etemp);

  // grab the data for this electron temperature region:
  std::vector<std::vector<std::vector<std::vector<double>>>> csk_data =
      Cdata->get_csk_data(i);

  // get the electron temperature eval points, too:
  std::vector<double> etemp_pts = Cdata->get_etemp_data(i);

  // make a vector for the return value:
  std::vector<std::vector<std::vector<double>>> interp_data(
      csk_data[0].size(),
      std::vector<std::vector<double>>(
          csk_data[0][0].size(),
          std::vector<double>(csk_data[0][0][0].size(), 0.0)));

  // we assume the matrix of csk data is not jagged...
  for (size_t a = 0; a < csk_data[0].size(); a++) {
    for (size_t b = 0; b < csk_data[0][0].size(); b++) {
      for (size_t c = 0; c < csk_data[0][0][0].size(); c++) {
        //stride through and get the correct "stripe" of data
        std::vector<double> interp_pts(csk_data.size(), 0.0);
        for (size_t d = 0; d < csk_data.size(); d++) {
          interp_pts[d] = csk_data[d][a][b][c];
        }
        //interpolate!
        interp_data[a][b][c] = interpolate_etemp(interp_pts, etemp_pts, etemp);
      }
    }
  }

  return interp_data;
}

// Use Lagrange interpolation on a set of csk_data points for some electron
// temperature
double ComptonInterp::interpolate_etemp(const std::vector<double> &csk_data,
                                        const std::vector<double> &etemp_data,
                                        const double etemp) {
  // Here, we check to see if the data point lies directly on one of our
  // eval points. If it does, we just return the CSK at the data point.

  double phi = 1.0;
  for (size_t j = 0; j < netemp_local; ++j) {
    if ((etemp == etemp_data[j])) {
      return csk_data[j];
    } else {
      phi *= (etemp - etemp_data[j]);
    }
  }
  double val = 0.;
  for (size_t j = 0; j < netemp_local; ++j) {
    val += csk_data[j] * (phi * prod_etemp[j]) / (etemp - etemp_data[j]);
  }

  return val;
}
}
