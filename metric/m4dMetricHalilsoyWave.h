//////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2018, USTC, Liang Yien
//
//  Einstein's Camera - A Relativity Theory Based Physical Renderer
//
//  A Program for the (Undergraduate Course) Optics Presentation
//
//  Modified from GeoVis  
//  And also contained Motion4D, TinyScheme
//  Uses GNU Scientific Library
//  Liang Yien, 
//  Student,
//  University of Science and Technology of China
// 
//  I was confused about the license of Motion4D and GeoVis, in their
//  code they claim they obey the GPL, while the author claim the 
//  software uses CPC License. Anyway, I would be safe to modified
//  the code for presentation in the class. But for the convience
//  of my classmate I might offer them my code , which is done rather
//  privately, like send them code in E-mail.I would immediately stop
//  if I am notified that it violates the author's willing.
//  
//  This code will be kept personal and not distribute publicly given
//  the fact that Motion4D and GeoVis were claimed to applied CPC
//  License. However, I would be gald to have my code in publid freely
//  using GPL if I am convicted that Motion4D and GeoVis is under GPL. 
//
//  Einstein's Camera was  written for my personally purpose only at 
//  the beginning and still is now WITHOUT ANY WARRANTY.No MERCHANTABILITY 
//  or FITNESS FOR A PARTICULAR PURPOSE is guaranteed. 
///////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2013, Universitaet Stuttgart, VISUS, Thomas Mueller
//
//  GeoViS is free software: you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  GeoViS is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GeoViS.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave
//
//  Motion4D- A library for numerical relativity calulation
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or (at
//  your option) any later version.
//  
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
////////////////////////////////////////////////////////////////////////////////////
// GSL - GNU SCIENTIFIC LIBRARY
// 
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Mark Galassi, James Theiler,
// Brian Gough, Gerard Jungman and many others
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
/////////////////////////////////////////////////////////////////////////////////////
// T I N Y S C H E M E    1 . 3 9
//   Dimitrios Souflis (dsouflis@acm.org)
//   Based on MiniScheme (original credits follow)
// (MINISCM)               coded by Atsushi Moriwaki (11/5/1989)
// (MINISCM)           E-MAIL :  moriwaki@kurims.kurims.kyoto-u.ac.jp
// (MINISCM) This version has been modified by R.C. Secrist.
// (MINISCM)
// (MINISCM) Mini-Scheme is now maintained by Akira KIDA.
// (MINISCM)
// (MINISCM) This is a revised and modified version by Akira KIDA.
// (MINISCM)	current version is 0.85k4 (15 May 1994)
/////////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------------
/*
    m4dMetricHalilsoyWave.h
 
  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave
 
 
    This file is part of the m4d-library.
 
    The m4d-library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    The m4d-library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.
 
 */

/*!  \class  m4d::MetricHalilsoyWave
     \brief  Gravitational wave solution given by Halilsoy.

       The Halilsoy standing gravitational wavein cylindrical coordinates (t,rho,phi,z) is given by
       \f[ ds^2 = e^{-2U}\left[e^{2k}\left(d\rho^2-dt^2\right) + \rho^2d\varphi^2\right] + e^{2U}\left(dz+A\,d\varphi\right)^2,\f]
       where \f$e^{-2U}=\cosh^2\alpha\,e^{-2CJ_0(\rho)\cos t}+\sinh^2\alpha\,e^{2CJ_0(\rho)\cos t}\f$, \f$A=-2C\sinh(2\alpha)\rho J_1(\rho)\sin t\f$,
       \f$k=\frac{1}{2}C^2\left[\rho^2\left(J_o(\rho)^2+J_1(\rho)^2\right)-2\rho J_0(\rho)J_1(\rho)\cos^2t\right]\f$, and the Bessel functions J_i.

     For further information see<br>
     M. Halilsoy, <b>Cross-Polarized Cylindrical Gravitational Waves of Einstein and Rosen</b>,<br>
     Il Nuovo Cimento <b>102 B</b>, 563 (1988).<br>
     or<br>
     Exact solutions (p.353).
 */
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_HALILSOY_WAVE_H
#define M4D_METRIC_HALILSOY_WAVE_H

#include "m4dMetric.h"
#include "gsl/gsl_sf_bessel.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   MetricHalilsoyWave
// ---------------------------------------------------
class MetricHalilsoyWave: public Metric
{
	public:
    MetricHalilsoyWave ( double alpha = 0.0, double C = 1.0 );
    ~MetricHalilsoyWave ();


// --------- public methods --------------
 public:
  virtual bool   calculateMetric        ( const double* pos );
  virtual bool   calculateChristoffels  ( const double* pos );
  virtual bool   calculateChrisD        ( const double* pos );

  virtual void   localToCoord           ( const double* pos, const double* ldir, double* dir,
                                          enum_nat_tetrad_type type =enum_nat_tetrad_default );
  virtual void   coordToLocal           ( const double* pos, const double* cdir, double* ldir,
                                          enum_nat_tetrad_type type = enum_nat_tetrad_default );

  virtual bool   breakCondition         ( const double* pos);

  virtual double testConstraint         ( const double y[], const double kappa );

  virtual bool   setParam               ( std::string pName, double val );

  virtual bool   report    ( const vec4 pos, const vec4 cdir, std::string &text );

// --------- protected methods -----------
 protected:
  virtual void setStandardValues ();

// ----- specific protected methods ------
  void  calcMetricFunc     ( const double* pos );
  void  calcDiffMetricFunc ( const double* pos );

// -------- protected attribute ----------
  protected:
    double mAlpha;
    double mC;
    
    // metric functions
    double fV, fA, fK;
    
    // derivatives of the metric functions
    double fV_t, fV_tt, fV_rho, fV_rhorho, fV_trho;
    double fA_t, fA_tt, fA_rho, fA_rhorho, fA_trho;
    double fK_t, fK_tt, fK_rho, fK_rhorho, fK_trho;
};

} // end namespace m4d

#endif
