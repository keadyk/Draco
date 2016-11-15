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

#include <vector>
#include <cstdlib>
#include <iostream>

namespace rtt_compton {

class ComptonData {
  private:
    //! Interpolation breakpoints for compton library:
    std::vector<double> etemp_breakpts;
    std::vector<double> gin_breakpts;
    std::vector<double> gout_breakpts;
    //! Evaluation points for the compton library:
    std::vector<double> etemp_pts;
    std::vector<double> gin_pts;
    std::vector<double> gout_pts;
    std::vector<double> xi_pts;

    //! flags for Lagrange interpolation and Legendre moment angular treatment
    bool lagrange, legendre;

    //! Raw CSK data in electron temperature, frequency in, frequency out,
    //! and xi (or Legendre moments of xi)
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

    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    // ACCESSORS -------------------------------------------------------------//
    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    //! Accessors for various data arrays
    size_t get_n_etemp_breakpts() { return etemp_breakpts.size(); }
    size_t get_n_gin_breakpts() { return gin_breakpts.size(); }
    size_t get_n_gout_breakpts() { return gout_breakpts.size(); }

    size_t get_n_etemp_pts() { return etemp_pts.size(); }
    size_t get_n_gin_pts() { return gin_pts.size(); }
    size_t get_n_gout_pts() { return gout_pts.size(); }
    size_t get_n_xi_pts() { return xi_pts.size(); }

    std::vector<double> get_etemp_breakpts() { return etemp_breakpts; }
    std::vector<double> get_gin_breakpts() { return gin_breakpts; }
    std::vector<double> get_gout_breakpts() { return gout_breakpts; }

    std::vector<double> get_etemp_pts() { return etemp_pts; }
    std::vector<double> get_gin_pts() { return gin_pts; }
    std::vector<double> get_gout_pts() { return gout_pts; }
    std::vector<double> get_xi_pts() { return xi_pts; }

    // (get ALL CSK data in one fell swoop)
    std::vector<std::vector<std::vector<std::vector<double>>>> get_csk_data()
                                                            { return csk_data; }

    //! Const accessors for various data arrays
    size_t get_n_etemp_breakpts() const { return etemp_breakpts.size(); }
    size_t get_n_gin_breakpts() const { return gin_breakpts.size(); }
    size_t get_n_gout_breakpts() const { return gout_breakpts.size(); }

    size_t get_n_etemp_pts() const { return etemp_pts.size(); }
    size_t get_n_gin_pts() const { return gin_pts.size(); }
    size_t get_n_gout_pts() const { return gout_pts.size(); }
    size_t get_n_xi_pts() const { return xi_pts.size(); }

    std::vector<double> get_etemp_breakpts() const { return etemp_breakpts; }
    std::vector<double> get_gin_breakpts() const { return gin_breakpts; }
    std::vector<double> get_gout_breakpts() const { return gout_breakpts; }

    std::vector<double> get_etemp_pts() const { return etemp_pts; }
    std::vector<double> get_gin_pts() const { return gin_pts; }
    std::vector<double> get_gout_pts() const { return gout_pts; }
    std::vector<double> get_xi_pts() const { return xi_pts; }

    // (get ALL CSK data in one fell swoop)
    std::vector<std::vector<std::vector<std::vector<double>>>> get_csk_data() 
                                                      const { return csk_data; }
    

    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    // SETTERS ---------------------------------------------------------------//
    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    
    void set_csk_data(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& data)
      { csk_data = data; }

    // set all evaluation points at once:
    void set_evalpts(const std::vector<double>& etempdata, 
                      const std::vector<double>& gindata, 
                      const std::vector<double>& goutdata, 
                      const std::vector<double>& xidata)
    {
      etemp_pts = etempdata;
      gin_pts = gindata;
      gout_pts = goutdata;
      xi_pts = xidata;
    }

    // set all break points at once:
    void set_breakpts(const std::vector<double>& etempdata, 
                      const std::vector<double>& gindata, 
                      const std::vector<double>& goutdata)
    {
      etemp_breakpts = etempdata;
      gin_breakpts = gindata;
      gout_breakpts = goutdata;
    }


    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    // BOOL ACCESSORS --------------------------------------------------------//
    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    //! retrieve bools:  
    bool is_lagrange() { return lagrange; }
    bool is_legendre() { return legendre; }

    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    // DATA SUBSET ACCESSORS -------------------------------------------------//
    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    // get all data for a particular etemp
    std::vector<double> get_etemp_data(const size_t region) const
    {
      size_t ppr_etemp = etemp_pts.size()/(etemp_breakpts.size()-1);
      std::vector<double>etemp_data(ppr_etemp, 0.0);
      size_t index1 = region * ppr_etemp;
      for(size_t a=0; a<ppr_etemp; a++)
      { etemp_data[a] = etemp_pts[index1 + a]; }
      return etemp_data;  
    }

    std::vector<double> get_etemp_data(const size_t region)
    {
      size_t ppr_etemp = etemp_pts.size()/(etemp_breakpts.size()-1);
      std::vector<double>etemp_data(ppr_etemp, 0.0);
      size_t index1 = region * ppr_etemp;
      for(size_t a=0; a<ppr_etemp; a++)
      { etemp_data[a] = etemp_pts[index1 + a]; }
      return etemp_data;  
    }

    // get all data for a particlar etemp, interpolation region, and x/y region
    std::vector<std::vector<std::vector<std::vector<double>>>>get_csk_data(
                                                               size_t etemp_reg,
                                                               size_t interp_reg,
                                                               size_t x_reg,
                                                               size_t y_reg) const
    {
      // use the provided regions to determine the proper indices:
      size_t ppr_etemp = etemp_pts.size()/(etemp_breakpts.size()-1);
      size_t ppr_gin = gin_pts.size()/(gin_breakpts.size()-1);
      size_t ppr_gout = gout_pts.size()/(gout_breakpts.size()-1);

      // use the interp region (1,2, or 3) to calculate the offset into the gin 
      // index:
      size_t offset = interp_reg*gin_pts.size();

      // allocate a return vector
      std::vector<std::vector<std::vector<std::vector<double>>>> select_data
      (ppr_etemp, std::vector<std::vector<std::vector<double>>>(ppr_gin, 
       std::vector<std::vector<double>>(ppr_gout, 
       std::vector<double>(xi_pts.size(), 0.0) ) ) );
      
      // fill it with data:
      for(size_t a=0; a<ppr_etemp; a++)
      {
        for(size_t b=0; b<ppr_gin; b++)
        {
          for(size_t c=0; c<ppr_gout; c++)
          {
            // fill in the whole xi row from the csk_data;
            select_data[a][b][c] = csk_data[ppr_etemp*etemp_reg + a]
                                           [offset + ppr_gin*x_reg + b]
                                           [ppr_gout*y_reg + c];
          } 
        }
      }
      return select_data;
    }

    // get all data for a particlar etemp, interpolation region, and x/y region
    std::vector<std::vector<std::vector<std::vector<double>>>>get_csk_data(
                                                               size_t etemp_reg,
                                                               size_t interp_reg,
                                                               size_t x_reg,
                                                               size_t y_reg)
    {
      // use the provided regions to determine the proper indices:
      size_t ppr_etemp = etemp_pts.size()/(etemp_breakpts.size()-1);
      size_t ppr_gin = gin_pts.size()/(gin_breakpts.size()-1);
      size_t ppr_gout = gout_pts.size()/(gout_breakpts.size()-1);

      // use the interp region (1,2, or 3) to calculate the offset into the gin 
      // index:
      size_t offset = interp_reg*gin_pts.size();

      // allocate a return vector
      std::vector<std::vector<std::vector<std::vector<double>>>> select_data
      (ppr_etemp, std::vector<std::vector<std::vector<double>>>(ppr_gin, 
       std::vector<std::vector<double>>(ppr_gout, 
       std::vector<double>(xi_pts.size(), 0.0) ) ) );
      
      // fill it with data:
      for(size_t a=0; a<ppr_etemp; a++)
      {
        for(size_t b=0; b<ppr_gin; b++)
        {
          for(size_t c=0; c<ppr_gout; c++)
          {
            // fill in the whole xi row from the csk_data;
            select_data[a][b][c] = csk_data[ppr_etemp*etemp_reg + a]
                                           [offset + ppr_gin*x_reg + b]
                                           [ppr_gout*y_reg + c];
          } 
        }
      }
      return select_data;
    }

    // get all data for a particlar etemp
    std::vector<std::vector<std::vector<std::vector<double>>>>get_csk_data(
                                                        size_t etemp_reg) const
    {
      // use the provided regions to determine the proper indices:
      size_t ppr_etemp = etemp_pts.size()/(etemp_breakpts.size()-1);
      // allocate a return vector
      std::vector<std::vector<std::vector<std::vector<double>>>> select_data
      (ppr_etemp, std::vector<std::vector<std::vector<double>>>(3*gin_pts.size(), 
       std::vector<std::vector<double>>(gout_pts.size(), 
       std::vector<double>(xi_pts.size(), 0.0) ) ) );
      
      // for each local etemp interpolation point in this breakpt region:
      for(size_t a=0; a<ppr_etemp; a++)
      {
        // fill in all x/y/xi data from the csk_data;
        select_data[a] = csk_data[ppr_etemp*etemp_reg + a];
      }
      return select_data;
    }

    // get all data for a particlar etemp
    std::vector<std::vector<std::vector<std::vector<double>>>>get_csk_data(
                                                              size_t etemp_reg)
    {
      // use the provided regions to determine the proper indices:
      size_t ppr_etemp = etemp_pts.size()/(etemp_breakpts.size()-1);
      // allocate a return vector
      std::vector<std::vector<std::vector<std::vector<double>>>> select_data
      (ppr_etemp, std::vector<std::vector<std::vector<double>>>(3*gin_pts.size(), 
       std::vector<std::vector<double>>(gout_pts.size(), 
       std::vector<double>(xi_pts.size(), 0.0) ) ) );
      
      // for each local etemp interpolation point in this breakpt region:
      for(size_t a=0; a<ppr_etemp; a++)
      {
        // fill in all x/y/xi data from the csk_data;
        select_data[a] = csk_data[ppr_etemp*etemp_reg + a];
      }
      return select_data;
    }

}; // end class compton_data

} // end namespace rtt_compton

#endif
