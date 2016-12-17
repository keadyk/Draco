//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/ComptonData.hh
 * \author Kendra Keady
 * \date   Tues Nov 8 2016 (Election Day!)
 * \brief  Header/Implementation file for ComptonData class (holds data for 
 *         Compton scattering interpolation)
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
#ifndef __compton_ComptonData_hh__
#define __compton_ComptonData_hh__

#include <cstdlib>
#include <iostream>
#include <vector>

namespace rtt_compton {

class ComptonData {
private:
  //! Interpolation breakpoints for multigroup compton library:
  std::vector<double> etemp_breakpts;
  std::vector<double> etemp_pts;

  //! Group boundaries of library:
  std::vector<double> group_bounds;

  //! number of legendre moments, frequency groups, interp points per etemp:
  size_t n_leg, n_grp;

  //! Raw CSK data in electron temperature, g in, g out,
  //! and Legendre moments of xi
  std::vector<std::vector<std::vector<std::vector<double>>>> csk_data;

public:
  //! Constructor
  ComptonData(const size_t netemp, const size_t ngrp, const size_t nleg)
      : n_leg(nleg), n_grp(ngrp) {
    // reserve memory for the data vectors whose sizes we now know:
    etemp_pts.reserve(netemp);
    group_bounds.reserve(ngrp + 1);
    // initialize the raw csk data vector
    csk_data.assign(netemp,
                    std::vector<std::vector<std::vector<double>>>(
                        ngrp, std::vector<std::vector<double>>(
                                  ngrp, std::vector<double>(nleg, 0.0))));
  }

  //! Destructor
  ~ComptonData() {}

  //------------------------------------------------------------------------//
  //------------------------------------------------------------------------//
  // ACCESSORS -------------------------------------------------------------//
  //------------------------------------------------------------------------//
  //------------------------------------------------------------------------//
  //! Accessors for various data arrays
  size_t get_n_etemp_breakpts() { return etemp_breakpts.size(); }

  size_t get_n_etemp_pts() { return etemp_pts.size(); }

  size_t get_n_leg() { return n_leg; }

  std::vector<double> get_etemp_breakpts() { return etemp_breakpts; }

  std::vector<double> get_etemp_pts() { return etemp_pts; }

  std::vector<double> get_group_bds() { return group_bounds; }

  // (get ALL CSK data in one fell swoop)
  std::vector<std::vector<std::vector<std::vector<double>>>> get_csk_data() {
    return csk_data;
  }

  //! Const accessors for various data arrays
  size_t get_n_etemp_breakpts() const { return etemp_breakpts.size(); }

  size_t get_n_etemp_pts() const { return etemp_pts.size(); }

  size_t get_n_leg() const { return n_leg; }

  std::vector<double> get_etemp_breakpts() const { return etemp_breakpts; }

  std::vector<double> get_etemp_pts() const { return etemp_pts; }

  std::vector<double> get_group_bds() const { return group_bounds; }

  // (get ALL CSK data in one fell swoop)
  std::vector<std::vector<std::vector<std::vector<double>>>>
  get_csk_data() const {
    return csk_data;
  }

  //------------------------------------------------------------------------//
  //------------------------------------------------------------------------//
  // SETTERS ---------------------------------------------------------------//
  //------------------------------------------------------------------------//
  //------------------------------------------------------------------------//

  void set_csk_data(
      const std::vector<std::vector<std::vector<std::vector<double>>>> &data) {
    csk_data = data;
  }

  void set_evalpts(const std::vector<double> &etempdata,
                   const std::vector<double> &groupdata) {
    etemp_pts = etempdata;
    group_bounds = groupdata;
  }

  void set_breakpts(const std::vector<double> &etempdata) {
    etemp_breakpts = etempdata;
  }

  //------------------------------------------------------------------------//
  //------------------------------------------------------------------------//
  // DATA SUBSET ACCESSORS -------------------------------------------------//
  //------------------------------------------------------------------------//
  //------------------------------------------------------------------------//
  // get all csk data for a particular etemp interpolation regions:
  std::vector<std::vector<std::vector<std::vector<double>>>>
  get_csk_data(const size_t region) const {
    // get the number of local electron temperature points:
    size_t ppr_etemp = etemp_pts.size() / (etemp_breakpts.size() - 1);
    std::vector<std::vector<std::vector<double>>> dummy;
    std::vector<std::vector<std::vector<std::vector<double>>>> csk_subset(
        ppr_etemp, dummy);
    size_t index1 = region * ppr_etemp;
    for (size_t a = 0; a < ppr_etemp; a++) {
      csk_subset[a] = csk_data[index1 + a];
    }
    return csk_subset;
  }

  std::vector<std::vector<std::vector<std::vector<double>>>>
  get_csk_data(const size_t region) {
    // get the number of local electron temperature points:
    size_t ppr_etemp = etemp_pts.size() / (etemp_breakpts.size() - 1);
    std::vector<std::vector<std::vector<double>>> dummy;
    std::vector<std::vector<std::vector<std::vector<double>>>> csk_subset(
        ppr_etemp, dummy);
    size_t index1 = region * ppr_etemp;
    for (size_t a = 0; a < ppr_etemp; a++) {
      csk_subset[a] = csk_data[index1 + a];
    }
    return csk_subset;
  }

  // get all etemp data for a particular interpolation regions:
  std::vector<double> get_etemp_data(const size_t region) const {
    size_t ppr_etemp = etemp_pts.size() / (etemp_breakpts.size() - 1);
    std::vector<double> etemp_data(ppr_etemp, 0.0);
    size_t index1 = region * ppr_etemp;
    for (size_t a = 0; a < ppr_etemp; a++) {
      etemp_data[a] = etemp_pts[index1 + a];
    }
    return etemp_data;
  }

  std::vector<double> get_etemp_data(const size_t region) {
    size_t ppr_etemp = etemp_pts.size() / (etemp_breakpts.size() - 1);
    std::vector<double> etemp_data(ppr_etemp, 0.0);
    size_t index1 = region * ppr_etemp;
    for (size_t a = 0; a < ppr_etemp; a++) {
      etemp_data[a] = etemp_pts[index1 + a];
    }
    return etemp_data;
  }

}; // end class compton_data

} // end namespace rtt_compton

#endif
