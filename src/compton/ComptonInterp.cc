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

namespace rtt_compton {

ComptonInterp::ComptonInterp(rtt_dsxx::SP<const ComptonData> Cdata_) :
                                                                Cdata(Cdata_)
{}
/*
    xs = xs0;	
    fs = fs0;	
    cs.resize(xs.size());
    n = xs.size();
	// Pre-compute cxs_{j} = product_{l \neq j} 1/(xs[j] - xs[l]) 
    for (int j = 0; j < n; ++j)
	{
		cs[j] = 1.;
	    for (int l = 0; l < n; ++l)
		{
			if (l != j)
		    {
			    cs[j] = cs[j]/(xs[j] - xs[l]);
			}
		}
    }
}
    */


/*
double ComptonInterp::interpolate_etemp(const double etemp)
{
  double i = find_local_region(etemp); 
	int n = xs.size();
    double phi = 1.;
    for (int j=0; j < n; ++j)
    {
        phi = (x-xs[j]) * phi;
    }
    double val = 0.;
    for (int j=0; j < n; ++j)
    {
        val = val + fs[j] * (phi*cs[j])/(x-xs[j]);
    }  
    return val;

}

double ComptonInterp::interpolate_gamma(const double gin, const double gout,
                                        const double xi)
{
  Region this_reg = find_global_region(gin, gout);
  
  std::pair<double, double>ij = find_local_region(x, y); 
}
*/

}
