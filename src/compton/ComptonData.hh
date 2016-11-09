//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/ComptonInterp.hh
 * \author Kendra Keady
 * \date   Tues Nov 8 2016 (Election Day!)
 * \brief  Header file for ComptonData class
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
#ifndef __compton_ComptonData_hh__
#define __compton_ComptonData_hh__

#include <vector>
#include <cstdlib>

namespace rtt_compton {

class ComptonData {
  private:
    //! Evaluation points for the compton library:
    std::vector<double> etemp_breakpts;
    std::vector<double> gin_breakpts;
    std::vector<double> gout_breakpts;
    std::vector<double> etemp_pts;
    std::vector<double> gin_pts;
    std::vector<double> gout_pts;
    std::vector<double> xi_pts;

    //! flags for lagrange interpolation, legendre moment treatment
    bool lagrange, legendre;

    // Raw CSK data read in from library file:
    std::vector<std::vector<std::vector<std::vector<double>>>> csk_data;

  public:
    //! Constructor
    ComptonData(const size_t netemp, const size_t ngin, const size_t ngout, 
                const size_t nxi, const bool lag=true, const bool leg=false) : 
                                                                  lagrange(lag),
                                                                  legendre(leg)
    {
      // reserve memory for the data vectors whose sizes we now know:
      etemp_pts.reserve(netemp);
      gin_pts.reserve(ngin);
      gout_pts.reserve(ngout);
      xi_pts.reserve(nxi);
      // initialize the raw csk data vector
      size_t dim = (lagrange) ? 3*ngin : ngin;
      csk_data.assign(netemp, 
                      std::vector<std::vector<std::vector<double>>>(dim, 
                      std::vector<std::vector<double>>(ngout, 
                      std::vector<double>(nxi, 0.0))));

    }

    //! Destructor
    ~ComptonData() {}

    //! Accessors for various data arrays
    std::vector<double> get_etemp_breakpts() { return etemp_breakpts; }
    std::vector<double> get_gin_breakpts() { return gin_breakpts; }
    std::vector<double> get_gout_breakpts() { return gout_breakpts; }

    std::vector<double> get_etemp_pts() { return etemp_pts; }
    std::vector<double> get_gin_pts() { return gin_pts; }
    std::vector<double> get_gout_pts() { return gout_pts; }
    std::vector<double> get_xi_pts() { return xi_pts; }

    std::vector<std::vector<std::vector<std::vector<double>>>> get_csk_data()
                                                            { return csk_data; }

    //! Set various data arrays
    void set_etemp_breakpts(const std::vector<double>& data)
      { etemp_breakpts = data; }
    void set_gin_breakpts(const std::vector<double>& data)
      { gin_breakpts = data; }
    void set_gout_breakpts(const std::vector<double>& data)
      { gout_breakpts = data; }
    void set_etemp_pts(const std::vector<double>& data)
      { etemp_pts = data; }
    void set_gin_pts(const std::vector<double>& data)
      { gin_pts = data; }
    void set_gout_pts(const std::vector<double>& data)
      { gout_pts = data; }
    void set_xi_pts(const std::vector<double>& data)
      { xi_pts = data; }

    void set_csk_data(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& data)
      { csk_data = data; }



};

}

#endif
