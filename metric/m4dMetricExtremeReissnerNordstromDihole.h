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
    m4dMetricExtremeReissnerNordstromDihole.h 
 
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

/*!  \class  m4d::MetricExtremeReissnerNordstromDihole
     \brief  Extreme Reissner-Nordstrom Dihole in Cartesian coordinates (t,x,y,z).
  
             The line element is given by
  
             \f[ds^2 = -\frac{dt}{U^2}+U^2\left(dx^2+dy^2+dz^2\right),\f]
  
             where \f$U=1+\frac{M_1}{r_1}+\frac{M_2}{r_2}\f$, \f$r_1=\sqrt{x^2+y^2+(z-1)^2}\f$, and \f$r_2=\sqrt{x^2+y^2+(z+1)^2}\f$.
  
             The metric is taken from<br>
             S. Chandrasekhar,<br>
             <b>The two-centre problem in general relativity: the scattering of radiation by two extreme Reissner-Nordstrom black-holes,</b><br>
             Proc. R. Soc. Lond. A <b>421</b>, 227-258 (1989).
  
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_EXTREME_REISSNER_DIHOLE_H
#define M4D_METRIC_EXTREME_REISSNER_DIHOLE_H

#include "m4dMetric.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

namespace m4d
{
 
// ---------------------------------------------------
//    class definition:   MetricExtremeReissnerNordstromDihole
// ---------------------------------------------------
class MetricExtremeReissnerNordstromDihole : public Metric
{
 public:
  MetricExtremeReissnerNordstromDihole ( double mass1 = 1.0, double mass2 = 1.0 );
  virtual ~MetricExtremeReissnerNordstromDihole (); 
  
// --------- public methods -----------
 public:
  virtual bool   calculateMetric        ( const double* pos );
  virtual bool   calculateChristoffels  ( const double* pos );
  virtual bool   calculateChrisD        ( const double* pos );
  
  virtual void   localToCoord           ( const double* pos, const double* ldir, double* dir,
                                          enum_nat_tetrad_type  type = enum_nat_tetrad_default );
  virtual void   coordToLocal           ( const double* pos, const double* cdir, double* ldir,
                                          enum_nat_tetrad_type  type = enum_nat_tetrad_default );
 
  virtual bool   breakCondition         ( const double* pos);

  virtual double testConstraint         ( const double y[], const double kappa );
 
  virtual bool   setParam               ( std::string pName, double val );

  virtual bool   transToTwoPlusOne      ( vec4 p, vec4 &cp );

  virtual bool   report ( const vec4 pos, const vec4 cdir, std::string &text );

// --------- specific public methods ----------
 public:
  bool  calcEscapeVelocityAlongZ ( const vec4 pos, double &betaEscape, double m1, double m2 );
  bool  calcEscapeVelocityAlongX ( const vec4 pos, double &betaEscape, double MASS );
  bool  calcPointOfEquilibrium ( double &pointOfEquilibrium );
  bool  calcASingularityPoint ( double &singularityPoint );
  bool  calcVelocityToEquilibriumPoint ( const vec4 pos, double &betaEquilibrium, double pointOfEquilibrium, double mass1, double mass2 );
  bool  calcPeriodicAlongX ( const vec4 pos, double &periodicAlongX, double MASS );
  bool  calcVelocityAtOriginAlongX ( const vec4 pos, double &betaAtOriginAlongX, double MASS );
  bool  calcPhotOrbitsXY ( double &r_Nm, double &r_Np, double mquer );
  bool  calcKsiCritXY ( const vec4 pos, double &ksicritXY, double r_Np, double MASS );
  bool  calcRadiusIscoXY ( double &radiusIscoXY, double MASS );
  bool  calcRadiusIucoXY ( double &radiusIucoXY, double MASS );
  bool  calcBetaCoXY ( double &betaCoXY, double radiusCoXY );
  bool  calcRadiusForOtherCo ( double &radiusForOtherCo, double z );
  bool  calcRadiusForOtherCoForUnequalMasses ( double &radiusForOtherCo, double z, double M1, double M2 );
  bool  calcHeightOfOtherPhotonOrbits ( double &heightOfOtherPhotonOrbits, double M );
  bool  calcVelocityForOtherTimelikeCo ( const vec4 pos, double &betaTimelike, double M );
  bool  calcHeightOfOtherPhotonOrbitsForUnequalMasses ( double &heightOfOtherPhotonOrbits1, double &heightOfOtherPhotonOrbits2, double M1, double M2 );
  bool  calcVelocityForOtherTimelikeCoForUnequalMasses ( const vec4 pos, double &betaTimelike, double M1, double M2 );

// --------- protected methods -----------
 protected:
  virtual void   setStandardValues ( );
          void   calc_r            ( const double y[] );
          double calc_U            ( const double y[] );
          void   calc_dU           ( const double y[], double &Ux, double &Uy, double &Uz );
          void   calc_ddU          ( const double y[], double &Uxx, double &Uyy, double &Uzz,
                                                       double &Uxy, double &Uxz, double &Uyz );

// -------- protected attribute ---------
 protected:
  double  r1, r2;
  double  mM1;
  double  mM2;

  gsl_integration_workspace* w;
  gsl_function F;
};

// ------- for some numerical calculations --------
/*!           ( no class members ! )           */

struct inflectionPoints_params { double M; };
double inflectionPoints (double x, void *params);
double inflectionPoints_deriv (double x, void *params);
void   inflectionPoints_fdf (double x, void *params, double *y, double *dy);

struct integrandForPeriodicAlongX_params { double M; double x0; };
double integrandForPeriodicAlongX (double x, void *params);

struct heightOfOtherPhotonOrbits_params { double M; };
double heightOfOtherPhotonOrbits (double z, void *params);
double heightOfOtherPhotonOrbits_deriv (double z, void *params);
void heightOfOtherPhotonOrbits_fdf (double z, void *params, double *y, double *dy);

struct heightOfOtherPhotonOrbitsForUnequalMasses_params { double M1,M2; };
double heightOfOtherPhotonOrbitsForUnequalMasses (double z, void *params);

} // end namespace m4d

#endif
