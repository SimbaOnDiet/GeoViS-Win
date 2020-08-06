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
   m4dMetricPlaneGravWave.h
 
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

/*!  \class  m4d::MetricPlaneGravWave
     \brief  Metric for a plane gravitational wave with finite logitudinal and
             infinite transversal extension (sandwich wave).

             With \f$ u=ct-x \f$ the line element is given by
             \f[ ds^2 = -c^2dt^2 + dx^2 + p^2\left(u\right)dy^2 + q^2\left(u\right)dz^2. \f]
             The functions \f$ p\left(u\right) \f$ and \f$ q\left(u\right) \f$
             are represented through
             \f[ p\left(u\right)=\left\{\begin{array}{lr}
                                          p_0 = \mathrm{const.} \quad &-a>u  \\
                                          L\left(u\right)\mathrm{e}^{m\left(u\right)} \quad  &-a\leq u \leq 0  \\
                                          1-u \quad &0<a
                                        \end{array}\right. \f]
             \f[ q\left(u\right)=\left\{\begin{array}{lr}
                                          q_0 = \mathrm{const.} \quad &-a>u  \\
                                          L\left(u\right)\mathrm{e}^{-m\left(u\right)} \quad  &-a\leq u \leq 0 \\
                                          1-u \quad &0<a
                                        \end{array}\right., \f]
             where the parameter \f$ a \f$ characterize the longitudinal
             extension of the wave. The Functions \f$ L\left(u\right) \f$ and
             \f$ m\left(u\right) \f$ are given by
             \f[ L\left(u\right) = 1 - u + \frac{1}{a^3}u^3 + \frac{1}{2a^3}u^4 \f]
             \f[ m\left(u\right) = \pm2\sqrt{3} \int \sqrt{\frac{u^{2} + au}{2a^{3}u - 2au^{3} - u^{4} - 2a^{3}}} \, \mathrm{d} u.\f]

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c}\partial_t,\quad \mathbf{e}_{(1)}=\partial_x,\quad \mathbf{e}_{(2)}= \frac{1}{p\left(u\right)}\partial_y,\quad \mathbf{e}_{(3)} = \frac{1}{p\left(u\right)}\partial_z.\f]

             Detailed discussions about the metric for a plane gravitational wave can be found
             in the standard literature, \ref lit_rindler "Rindler".
 */
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_PLANE_GRAV_WAVE_H
#define M4D_METRIC_PLANE_GRAV_WAVE_H

#include "m4dMetric.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   MetricPlaneGraveWave
// ---------------------------------------------------
class MetricPlaneGravWave: public Metric
{
 public:
  MetricPlaneGravWave ( double longExt = 1.0, double degree = 30.0 );
  virtual ~MetricPlaneGravWave ();

// ----------- public methods ------------
 public:
  virtual bool   calculateMetric        ( const double* pos );
  virtual bool   calculateChristoffels  ( const double* pos );
  virtual bool   calculateChrisD        ( const double* pos );

  virtual void   localToCoord           ( const double* pos, const double* ldir, double* dir,
                                          enum_nat_tetrad_type type = enum_nat_tetrad_default );
  virtual void   coordToLocal           ( const double* pos, const double* cdir, double* ldir,
                                          enum_nat_tetrad_type type = enum_nat_tetrad_default );

  virtual bool   breakCondition         ( const double* pos );

  virtual double testConstraint         ( const double y[], const double kappa );

  virtual bool   setParam               ( std::string pName, double val );
  virtual bool   transToTwoPlusOne      ( vec4 p, vec4 & cp );

  virtual bool   report  ( const vec4 pos, const vec4 cdir, std::string &text );

// --------- protected methods -----------
 protected:
  virtual void setStandardValues	();

// ----- specific protected methods ------
  double getValP           ( const double* pos );
  double getValQ           ( const double* pos );
  double getValDP          ( const double* pos );
  double getValDQ          ( const double* pos );
  double getValDDP         ( const double* pos );
  double getValDDQ         ( const double* pos );
  void calcFourierCoeff	   ( );

// --------- protected attribute ----------
 protected:
  bool dataCalculated;
  double mDegree;
  double mLongExt;
  double intConst;
  double p0;
  double q0;
  std::vector<double> fCoeffB;

  double *uData, *dmData;
};

} // end namespace m4d

#endif
