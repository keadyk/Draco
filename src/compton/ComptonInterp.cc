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
#include <sstream>

namespace rtt_compton {

ComptonInterp::ComptonInterp(rtt_dsxx::SP<const ComptonData> Cdata_)
    : Cdata(Cdata_) {
  xs = Cdata->get_gin_pts();
  ys = Cdata->get_gout_pts();
  etemps = Cdata->get_etemp_pts();

  nx_local = Cdata->get_n_gin_pts() / (Cdata->get_n_gin_breakpts() - 1);
  ny_local = Cdata->get_n_gout_pts() / (Cdata->get_n_gout_breakpts() - 1);
  netemp_local = Cdata->get_n_etemp_pts() / (Cdata->get_n_etemp_breakpts() - 1);

  nx_break = Cdata->get_n_gin_breakpts();
  ny_break = Cdata->get_n_gout_breakpts();
  netemp_break = Cdata->get_n_etemp_breakpts();

  prod_x.assign(xs.size(), 0.0);
  prod_y.assign(ys.size(), 0.0);
  prod_etemp.assign(etemps.size(), 0.0);

  // Pre-compute cxs_{i,j} = product_{l \neq j} 1/(xs[i,j] - xs[i,l])
  for (size_t i = 0; i < nx_break - 1; ++i) {
    for (size_t j = 0; j < nx_local; ++j) {
      size_t k = j + nx_local * i;
      prod_x[k] = 1.;
      for (size_t l = 0; l < nx_local; ++l) {
        if (l != j) {
          size_t k1 = l + nx_local * i;
          prod_x[k] = prod_x[k] / (xs[k] - xs[k1]);
        }
      }
    }
  }
  // Pre-compute cys_{i,j} = product_{l \neq j} 1/(ys[i,j] - ys[i,l])
  for (size_t i = 0; i < ny_break - 1; ++i) {
    for (size_t j = 0; j < ny_local; ++j) {
      size_t k = j + ny_local * i;
      prod_y[k] = 1.;
      for (size_t l = 0; l < ny_local; ++l) {
        if (l != j) {
          size_t k1 = l + ny_local * i;
          prod_y[k] = prod_y[k] / (ys[k] - ys[k1]);
        }
      }
    }
  }
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

Region ComptonInterp::find_global_region(const double gin, const double gout) {

  if (gout < gin && gout > (gin / (2 * gin + 1))) {
    // middle region!!
    return Region::MID;
  } else if (gout >= gin) {
    // top region!
    return Region::TOP;
  } else if (gout <= (gin / (2 * gin + 1))) {
    // bottom region!
    return Region::BOTTOM;
  } else {
    std::cout << "Uh-oh! This gamma in/out pair is invalid! " << std::endl;
    return Region::NONE;
  }
}

size_t ComptonInterp::find_etemp_region(const double etemp) {
  std::vector<double> breakpt_data = Cdata->get_etemp_breakpts();
  size_t iregion = binary_search(etemp, breakpt_data);

  return iregion;
}

void ComptonInterp::set_xy_and_region(const double gin, const double gout,
                                      const Region region,
                                      std::pair<size_t, size_t> &ij,
                                      std::pair<double, double> &xy) {

  // get x/y interpolation breakpoints from compton_data:
  std::vector<double> xs = Cdata->get_gin_breakpts();
  std::vector<double> ys = Cdata->get_gout_breakpts();

  // transformed gin/gout values:
  double x, y;
  // x is always equal to gin, so we can complete this binary search first:
  x = gin;
  xy.first = x;

  // use x value to determine breakpoint index:
  ij.first = binary_search(x, xs);

  switch (region) {
  case Region::BOTTOM:
    y = 1 - (2 * gin * gout) / (gin - gout);
    break;
  case Region::MID:
    y = 1 - (gin - gout) / (2 * gin * gout);
    break;
  case Region::TOP:
    y = gout - gin;
    break;
  default:
    std::ostringstream message;
    message << "ComptonInterp: gin/gout pair outside interpolation range";
    throw std::out_of_range(message.str());
    break;
  }
  // use y value to determine breakpoint index:
  xy.second = y;
  ij.second = binary_search(y, ys);
}

// binary search for grid index
// TODO: could be templated if we need more than just the double version
size_t ComptonInterp::binary_search(const double value,
                                    const std::vector<double> &bp_data) {
  // check that the value is actually between the bounds:
  if (value < bp_data[0] || value > bp_data.back()) {
    std::ostringstream message;
    message << "ComptonInterp: value outside binary search bounds";
    throw std::out_of_range(message.str());
  }

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
  while (abs(high - low) > 1) {
    index = bp_data.size() / 2;
    if (value < bp_data[index]) {
      high = index;
    } else {
      low = index;
    }
  }
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

  // for all of the x points (x3 regions)
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

// interpolate for a gin/gout pair, given all (possibly pre-interpolated)
// CSK data at the current electron temperature
// TODO: the inner interpolate_gin_gout call recalculates the region (and x/y),
// which is wasteful since gin/gout are fixed (and thus x/y and the region are
// also fixed)
std::vector<double> ComptonInterp::interpolate_gin_gout(
    const double gin, const double gout,
    const std::vector<std::vector<std::vector<double>>> &csk_data) {
  // get the electron temperature eval points, too:
  std::vector<double> x_pts = Cdata->get_gin_pts();
  std::vector<double> y_pts = Cdata->get_gout_pts();

  // make a vector for the return values:
  std::vector<double> interp_data(csk_data[0][0].size(), 0.0);

  // for each xi point
  for (size_t a = 0; a < interp_data.size(); a++) {
    //stride through and get the correct "stripe" of data
    std::vector<std::vector<double>> interp_pts(
        csk_data.size(), std::vector<double>(csk_data[0].size(), 0.0));
    for (size_t b = 0; b < csk_data.size(); b++) {
      for (size_t c = 0; c < csk_data[0].size(); c++) {
        interp_pts[b][c] = csk_data[b][c][a];
      }
    }
    //interpolate!
    interp_data[a] = interpolate_gin_gout(interp_pts, x_pts, y_pts, gin, gout);
  }

  return interp_data;
}

// Use Lagrange interpolation on a set of csk_data points for some electron
// temperature
double ComptonInterp::interpolate_etemp(const std::vector<double> &csk_data,
                                        const std::vector<double> &etemp_data,
                                        const double etemp) {
  double phi = 1.0;
  for (size_t j = 0; j < netemp_local; ++j) {
    phi *= (etemp - etemp_data[j]);
  }
  double val = 0.;
  for (size_t j = 0; j < netemp_local; ++j) {
    val += csk_data[j] * (phi * prod_etemp[j]) / (etemp - etemp_data[j]);
  }
  return val;
}

// Use 2-D Lagrange interpolation onf a set of csk_data points, for some gamma
// in/gamma out combination
double ComptonInterp::interpolate_gin_gout(
    const std::vector<std::vector<double>> &csk_data,
    const std::vector<double> &x_data, const std::vector<double> &y_data,
    const double gin, const double gout) {
  // First, we use gamma in and out to determine the global interpolation
  // region (bottom, middle, or top) bounded by the boundary-layer curves
  Region region = find_global_region(gin, gout);

  // find local intervals on which to do interpolation
  std::pair<size_t, size_t> ij;
  std::pair<double, double> xy;

  set_xy_and_region(gin, gout, region, ij, xy);
  size_t i_break = ij.first;
  size_t j_break = ij.second;
  double x = xy.first;
  double y = xy.second;

  // calculate offsets into the provided csk data, depending on the
  // global region (top, middle, or bottom, relative to boundary-layers in gin
  // and gout) and local region (breakpoint region in x/y)
  size_t x_offset = nx_local * i_break;
  size_t y_offset = ny_local * j_break;

  double phi1 = 1.0;
  double phi2 = 1.0;
  for (size_t j = 0; j < nx_local; ++j) {
    phi1 *= (x - x_data[x_offset + j]);
  }
  for (size_t j = 0; j < ny_local; ++j) {
    phi2 *= (y - y_data[y_offset + j]);
  }

  double val = 0.;
  for (size_t i_loc = 0; i_loc < nx_local; ++i_loc) {
    for (size_t j_loc = 0; j_loc < ny_local; ++j_loc) {
      size_t x_index = x_offset + i_loc;
      size_t y_index = y_offset + j_loc;

      double f_val = csk_data[int(region) * x_data.size() + x_index][y_index];
      val += f_val * (phi1 * prod_x[x_index]) * (phi2 * prod_y[y_index]) /
             ((x - x_data[x_index]) * (y - y_data[y_index]));
    }
  }

  return val;
}
}
