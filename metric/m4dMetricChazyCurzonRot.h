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
    m4dMetricChazyCurzonRot.h
 
  Copyright (c) 2011-2014  Thomas Mueller
 
 
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

/*!  \class  m4d::MetricChazyCurzonRot
     \brief  ChazyCurzonRot metric in Weyl coordinates (t,rho,phi,z).
  
             The line element is given by
  
             \f[ds^2 = e^{-2U}\left[e^{2k}\left(d\rho^2+dz^2\right)+\rho^2d\varphi^2\right]-e^{2U}\left(dt+A\,d\varphi\right)^2,\f]
             where \f$e^{-2U}=\cosh\frac{2m}{r}-p\sinh\frac{2m}{r}\f$, \f$2k=-m^2\rho^2/r^4\f$, \f$A=2qmz/r\f$, and \f$r^2=\rho^2+z^2\f$.
             The parameters \f$p,q\f$ are related via \f$p^2+q^2=1\f$.

             The natural local tetrad reads:


             See also
             Stephani, H., Kramer, D., MacCallum, M., Hoenselaers, C. and Herlt, E.,<br>
             Exact Solutions of the Einstein Field Equations<br>
             (Cambridge University Press, 2. edition, 2009)
  
*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_CHAZY_CURZON_ROT_H
#define M4D_METRIC_CHAZY_CURZON_ROT_H

#include "m4dMetric.h"

namespace m4d
{
 
// ---------------------------------------------------
//    class definition:   MetricChazyCurzonRot
// ---------------------------------------------------
class MetricChazyCurzonRot : public Metric
{
 public:
  MetricChazyCurzonRot ( double mass = 1.0, double p = 0.5 );
  virtual ~MetricChazyCurzonRot ();
  
// --------- public methods -----------
 public:
  virtual bool   calculateMetric        ( const double* pos );
  virtual bool   calculateChristoffels  ( const double* pos );
  virtual bool   calculateChrisD        ( const double* pos );
  
  virtual void   localToCoord           ( const double* pos, const double* ldir, double* dir, 
                                          enum_nat_tetrad_type  type = enum_nat_tetrad_cylinder );
  virtual void   coordToLocal           ( const double* pos, const double* cdir, double* ldir, 
                                          enum_nat_tetrad_type  type = enum_nat_tetrad_cylinder );
 
  virtual bool   breakCondition         ( const double* pos);

  virtual double testConstraint         ( const double y[], const double kappa );
 
  virtual bool   setParam ( std::string pName, double val );
  
  virtual bool   report ( const vec4 pos, const vec4 cdir, std::string &text );
  

// --------- protected methods -----------
 protected:
  virtual void setStandardValues ( );

          void calcUkA   ( const double* pos, double &em2U, double &k, double &A );
          void calcDUka  ( const double* pos, double &dUdrho, double &dUdz,
                                              double &dkdrho, double &dkdz,
                                              double &dAdrho, double &dAdz );
   
// -------- protected attribute ---------
 protected:
  double mMass;
  double m_p;
  double m_q;

};

} // end namespace m4d

#endif
