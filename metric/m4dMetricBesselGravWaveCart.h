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
   m4dMetricBesselGravWaveCart.h
 
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

/*!  \class  m4d::MetricBesselGravWaveCart
     \brief  Gravitational Bessel-Wave metric in cartesian coordinates (t,x,y,z).

             The line element is given by
             \f[ \mathrm{d} s^{2} = -\mathrm{e}^{2\left(K-U\right)} \mathrm{d}t^{2}
                                    + \frac{\mathrm{e}^{-2U}}{x^{2}+y^{2}} \left[
                                    \left(\mathrm{e}^{2K}x^{2}+y^{2}\right) \mathrm{d}x^{2}
                                    + 2xy\left(\mathrm{e}^{2K}-1\right) \mathrm{d}x \mathrm{d}y
                                    + \left(x^{2}+\mathrm{e}^{2K}y^{2}\right) \mathrm{d}y^{2} \right]
                                    + \mathrm{e}^{2U} \mathrm{d}z^{2}. \f]
             With \f$ \rho=\sqrt{x^2+y^2} \f$ the Functions \f$ U \f$ and \f$ K \f$
             are given by
             \f[ U = CJ_{0}\left(\rho\right)\cos\left(t\right) \f]
             \f[ K = \frac{1}{2}C^{2}\rho\left\{\rho\right[J_{0}\left(\rho\right)^{2}+J_{1}\left(\rho\right)^{2}\left]
                     -2J_{0}\left(\rho\right)J_{1}\left(\rho\right)\cos^{2}\left(t\right)\right\}, \f]
             where \f$ J_{n} \f$ are the Bessel functions of the first kind and
             \f$ C \f$ is a free Parameter.
             The metric is taken from
             D. Kramer, "An exact gravitational wave solution", Ann. Phys. Volume 9 (2000), 331-334,
             with the modification \f$ \left(\omega t-kz\right) \rightarrow t \f$,
             \f$ \left(kt-\omega z\right) \rightarrow z\f$ and a transformation
             to cartesian coordinates.
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_BESSELGRAVWAVECART_H
#define M4D_METRIC_BESSELGRAVWAVECART_H

#include "m4dMetric.h"
#include "gsl/gsl_sf_bessel.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   MetricBesselGravWaveCart
// ---------------------------------------------------
class MetricBesselGravWaveCart : public Metric
{
 public:
  MetricBesselGravWaveCart ( double C = 1.0 );
  virtual ~MetricBesselGravWaveCart ();
  
// --------- public methods -----------
 public:
  virtual bool   calculateMetric        ( const double* pos );
  virtual bool   calculateChristoffels  ( const double* pos );
  virtual bool	 calculateChrisD        ( const double* pos );
  
  virtual void   localToCoord           ( const double* pos, const double* ldir, double* dir,
                                          enum_nat_tetrad_type  type = enum_nat_tetrad_default );
  virtual void   coordToLocal           ( const double* pos, const double* cdir, double* ldir,
                                          enum_nat_tetrad_type  type = enum_nat_tetrad_default );
 
  virtual bool   breakCondition         ( const double* pos);

  virtual double testConstraint         ( const double y[], const double kappa );
 
  virtual bool   setParam               ( std::string pName, double val );

  virtual bool   transToTwoPlusOne      ( vec4 p, vec4 &cp );

  virtual bool   report       ( const vec4 pos, const vec4 cdir, std::string &text );

// --------- protected methods -----------
 protected:
  virtual void setStandardValues ();

// ----- specific protected methods ------
  void calcLTcoeffs ( const double* pos );
 
// -------- protected attribute ----------
 protected:
  double mC;
  // Coefficients for local Tetrad
  double ltA;
  double ltB;
  double ltC;
  double ltD;
  double ltE;
  double ltF;
  double ltT;
  double ltZ;

};

} // end namespace m4d

#endif
