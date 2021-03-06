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
   m4dMetricKasner.cpp
 
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
// -------------------------------------------------------------------------------

#include "m4dMetricKasner.h"

namespace m4d
{

#define eps 1.0e-6


/*! Standard constructor for the Kottler metric.
 *
 * \param  u : Khalatnikov-Lifshitz parameter.
 */
MetricKasner :: MetricKasner ( double u )
{
  mMetricName  = "Kasner";
  setCoordType(enum_coordinate_cartesian);
    
  mPhysicalUnits = enum_physical_constants_geom;
  mSpeedOfLight = 1.0;
  mGravConstant = 1.0;
  
  addParam("u",u); 
  mU = u;
  calc_parameters();

  setStandardValues();
}

MetricKasner :: ~MetricKasner()
{
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricKasner :: calculateMetric  ( const double* pos )
{
  double t = pos[0];

  double t1 = pow(t,p1);
  double t2 = t1*t1;
  double t3 = pow(t,p2);
  double t4 = t3*t3;
  double t5 = pow(t,p3);
  double t6 = t5*t5;

      g_compts[0][0] = -1.0;
      g_compts[0][1] = 0.0;
      g_compts[0][2] = 0.0;
      g_compts[0][3] = 0.0;
      g_compts[1][0] = 0.0;
      g_compts[1][1] = t2;
      g_compts[1][2] = 0.0;
      g_compts[1][3] = 0.0;
      g_compts[2][0] = 0.0;
      g_compts[2][1] = 0.0;
      g_compts[2][2] = t4;
      g_compts[2][3] = 0.0;
      g_compts[3][0] = 0.0;
      g_compts[3][1] = 0.0;
      g_compts[3][2] = 0.0;
      g_compts[3][3] = t6;

  return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricKasner :: calculateChristoffels  ( const double* pos )
{
  double t = pos[0];

  double t1 = 1/t;
  double t2 = p1*t1;
  double t3 = p2*t1;
  double t4 = p3*t1;
  double t5 = pow(t,p1);
  double t6 = t5*t5;
  double t9 = pow(t,p2);
  double t10 = t9*t9;
  double t13 = pow(t,p3);
  double t14 = t13*t13;

      christoffel[0][0][0] = 0.0;
      christoffel[0][0][1] = 0.0;
      christoffel[0][0][2] = 0.0;
      christoffel[0][0][3] = 0.0;
      christoffel[0][1][0] = 0.0;
      christoffel[0][1][1] = t2;
      christoffel[0][1][2] = 0.0;
      christoffel[0][1][3] = 0.0;
      christoffel[0][2][0] = 0.0;
      christoffel[0][2][1] = 0.0;
      christoffel[0][2][2] = t3;
      christoffel[0][2][3] = 0.0;
      christoffel[0][3][0] = 0.0;
      christoffel[0][3][1] = 0.0;
      christoffel[0][3][2] = 0.0;
      christoffel[0][3][3] = t4;
      christoffel[1][0][0] = 0.0;
      christoffel[1][0][1] = t2;
      christoffel[1][0][2] = 0.0;
      christoffel[1][0][3] = 0.0;
      christoffel[1][1][0] = t6*p1*t1;
      christoffel[1][1][1] = 0.0;
      christoffel[1][1][2] = 0.0;
      christoffel[1][1][3] = 0.0;
      christoffel[1][2][0] = 0.0;
      christoffel[1][2][1] = 0.0;
      christoffel[1][2][2] = 0.0;
      christoffel[1][2][3] = 0.0;
      christoffel[1][3][0] = 0.0;
      christoffel[1][3][1] = 0.0;
      christoffel[1][3][2] = 0.0;
      christoffel[1][3][3] = 0.0;
      christoffel[2][0][0] = 0.0;
      christoffel[2][0][1] = 0.0;
      christoffel[2][0][2] = t3;
      christoffel[2][0][3] = 0.0;
      christoffel[2][1][0] = 0.0;
      christoffel[2][1][1] = 0.0;
      christoffel[2][1][2] = 0.0;
      christoffel[2][1][3] = 0.0;
      christoffel[2][2][0] = t10*p2*t1;
      christoffel[2][2][1] = 0.0;
      christoffel[2][2][2] = 0.0;
      christoffel[2][2][3] = 0.0;
      christoffel[2][3][0] = 0.0;
      christoffel[2][3][1] = 0.0;
      christoffel[2][3][2] = 0.0;
      christoffel[2][3][3] = 0.0;
      christoffel[3][0][0] = 0.0;
      christoffel[3][0][1] = 0.0;
      christoffel[3][0][2] = 0.0;
      christoffel[3][0][3] = t4;
      christoffel[3][1][0] = 0.0;
      christoffel[3][1][1] = 0.0;
      christoffel[3][1][2] = 0.0;
      christoffel[3][1][3] = 0.0;
      christoffel[3][2][0] = 0.0;
      christoffel[3][2][1] = 0.0;
      christoffel[3][2][2] = 0.0;
      christoffel[3][2][3] = 0.0;
      christoffel[3][3][0] = t14*p3*t1;
      christoffel[3][3][1] = 0.0;
      christoffel[3][3][2] = 0.0;
      christoffel[3][3][3] = 0.0;

  return true;
}

/*! Calculate Jacobi matrix. 
 *
 *  \param pos : pointer to position.
 */
bool 
MetricKasner :: calculateChrisD ( const double* pos )
{
  double t = pos[0];

  double t1 = t*t;
  double t2 = 1/t1;
  double t3 = p1*t2;
  double t4 = p2*t2;
  double t5 = p3*t2;
  double t6 = pow(t,p1);
  double t7 = t6*t6;
  double t13 = pow(t,p2);
  double t14 = t13*t13;
  double t20 = pow(t,p3);
  double t21 = t20*t20;

      chrisD[0][0][0][0] = 0.0;
      chrisD[0][0][0][1] = 0.0;
      chrisD[0][0][0][2] = 0.0;
      chrisD[0][0][0][3] = 0.0;
      chrisD[0][0][1][0] = 0.0;
      chrisD[0][0][1][1] = 0.0;
      chrisD[0][0][1][2] = 0.0;
      chrisD[0][0][1][3] = 0.0;
      chrisD[0][0][2][0] = 0.0;
      chrisD[0][0][2][1] = 0.0;
      chrisD[0][0][2][2] = 0.0;
      chrisD[0][0][2][3] = 0.0;
      chrisD[0][0][3][0] = 0.0;
      chrisD[0][0][3][1] = 0.0;
      chrisD[0][0][3][2] = 0.0;
      chrisD[0][0][3][3] = 0.0;
      chrisD[0][1][0][0] = 0.0;
      chrisD[0][1][0][1] = 0.0;
      chrisD[0][1][0][2] = 0.0;
      chrisD[0][1][0][3] = 0.0;
      chrisD[0][1][1][0] = -t3;
      chrisD[0][1][1][1] = 0.0;
      chrisD[0][1][1][2] = 0.0;
      chrisD[0][1][1][3] = 0.0;
      chrisD[0][1][2][0] = 0.0;
      chrisD[0][1][2][1] = 0.0;
      chrisD[0][1][2][2] = 0.0;
      chrisD[0][1][2][3] = 0.0;
      chrisD[0][1][3][0] = 0.0;
      chrisD[0][1][3][1] = 0.0;
      chrisD[0][1][3][2] = 0.0;
      chrisD[0][1][3][3] = 0.0;
      chrisD[0][2][0][0] = 0.0;
      chrisD[0][2][0][1] = 0.0;
      chrisD[0][2][0][2] = 0.0;
      chrisD[0][2][0][3] = 0.0;
      chrisD[0][2][1][0] = 0.0;
      chrisD[0][2][1][1] = 0.0;
      chrisD[0][2][1][2] = 0.0;
      chrisD[0][2][1][3] = 0.0;
      chrisD[0][2][2][0] = -t4;
      chrisD[0][2][2][1] = 0.0;
      chrisD[0][2][2][2] = 0.0;
      chrisD[0][2][2][3] = 0.0;
      chrisD[0][2][3][0] = 0.0;
      chrisD[0][2][3][1] = 0.0;
      chrisD[0][2][3][2] = 0.0;
      chrisD[0][2][3][3] = 0.0;
      chrisD[0][3][0][0] = 0.0;
      chrisD[0][3][0][1] = 0.0;
      chrisD[0][3][0][2] = 0.0;
      chrisD[0][3][0][3] = 0.0;
      chrisD[0][3][1][0] = 0.0;
      chrisD[0][3][1][1] = 0.0;
      chrisD[0][3][1][2] = 0.0;
      chrisD[0][3][1][3] = 0.0;
      chrisD[0][3][2][0] = 0.0;
      chrisD[0][3][2][1] = 0.0;
      chrisD[0][3][2][2] = 0.0;
      chrisD[0][3][2][3] = 0.0;
      chrisD[0][3][3][0] = -t5;
      chrisD[0][3][3][1] = 0.0;
      chrisD[0][3][3][2] = 0.0;
      chrisD[0][3][3][3] = 0.0;
      chrisD[1][0][0][0] = 0.0;
      chrisD[1][0][0][1] = 0.0;
      chrisD[1][0][0][2] = 0.0;
      chrisD[1][0][0][3] = 0.0;
      chrisD[1][0][1][0] = -t3;
      chrisD[1][0][1][1] = 0.0;
      chrisD[1][0][1][2] = 0.0;
      chrisD[1][0][1][3] = 0.0;
      chrisD[1][0][2][0] = 0.0;
      chrisD[1][0][2][1] = 0.0;
      chrisD[1][0][2][2] = 0.0;
      chrisD[1][0][2][3] = 0.0;
      chrisD[1][0][3][0] = 0.0;
      chrisD[1][0][3][1] = 0.0;
      chrisD[1][0][3][2] = 0.0;
      chrisD[1][0][3][3] = 0.0;
      chrisD[1][1][0][0] = t7*p1*(2.0*p1-1.0)*t2;
      chrisD[1][1][0][1] = 0.0;
      chrisD[1][1][0][2] = 0.0;
      chrisD[1][1][0][3] = 0.0;
      chrisD[1][1][1][0] = 0.0;
      chrisD[1][1][1][1] = 0.0;
      chrisD[1][1][1][2] = 0.0;
      chrisD[1][1][1][3] = 0.0;
      chrisD[1][1][2][0] = 0.0;
      chrisD[1][1][2][1] = 0.0;
      chrisD[1][1][2][2] = 0.0;
      chrisD[1][1][2][3] = 0.0;
      chrisD[1][1][3][0] = 0.0;
      chrisD[1][1][3][1] = 0.0;
      chrisD[1][1][3][2] = 0.0;
      chrisD[1][1][3][3] = 0.0;
      chrisD[1][2][0][0] = 0.0;
      chrisD[1][2][0][1] = 0.0;
      chrisD[1][2][0][2] = 0.0;
      chrisD[1][2][0][3] = 0.0;
      chrisD[1][2][1][0] = 0.0;
      chrisD[1][2][1][1] = 0.0;
      chrisD[1][2][1][2] = 0.0;
      chrisD[1][2][1][3] = 0.0;
      chrisD[1][2][2][0] = 0.0;
      chrisD[1][2][2][1] = 0.0;
      chrisD[1][2][2][2] = 0.0;
      chrisD[1][2][2][3] = 0.0;
      chrisD[1][2][3][0] = 0.0;
      chrisD[1][2][3][1] = 0.0;
      chrisD[1][2][3][2] = 0.0;
      chrisD[1][2][3][3] = 0.0;
      chrisD[1][3][0][0] = 0.0;
      chrisD[1][3][0][1] = 0.0;
      chrisD[1][3][0][2] = 0.0;
      chrisD[1][3][0][3] = 0.0;
      chrisD[1][3][1][0] = 0.0;
      chrisD[1][3][1][1] = 0.0;
      chrisD[1][3][1][2] = 0.0;
      chrisD[1][3][1][3] = 0.0;
      chrisD[1][3][2][0] = 0.0;
      chrisD[1][3][2][1] = 0.0;
      chrisD[1][3][2][2] = 0.0;
      chrisD[1][3][2][3] = 0.0;
      chrisD[1][3][3][0] = 0.0;
      chrisD[1][3][3][1] = 0.0;
      chrisD[1][3][3][2] = 0.0;
      chrisD[1][3][3][3] = 0.0;
      chrisD[2][0][0][0] = 0.0;
      chrisD[2][0][0][1] = 0.0;
      chrisD[2][0][0][2] = 0.0;
      chrisD[2][0][0][3] = 0.0;
      chrisD[2][0][1][0] = 0.0;
      chrisD[2][0][1][1] = 0.0;
      chrisD[2][0][1][2] = 0.0;
      chrisD[2][0][1][3] = 0.0;
      chrisD[2][0][2][0] = -t4;
      chrisD[2][0][2][1] = 0.0;
      chrisD[2][0][2][2] = 0.0;
      chrisD[2][0][2][3] = 0.0;
      chrisD[2][0][3][0] = 0.0;
      chrisD[2][0][3][1] = 0.0;
      chrisD[2][0][3][2] = 0.0;
      chrisD[2][0][3][3] = 0.0;
      chrisD[2][1][0][0] = 0.0;
      chrisD[2][1][0][1] = 0.0;
      chrisD[2][1][0][2] = 0.0;
      chrisD[2][1][0][3] = 0.0;
      chrisD[2][1][1][0] = 0.0;
      chrisD[2][1][1][1] = 0.0;
      chrisD[2][1][1][2] = 0.0;
      chrisD[2][1][1][3] = 0.0;
      chrisD[2][1][2][0] = 0.0;
      chrisD[2][1][2][1] = 0.0;
      chrisD[2][1][2][2] = 0.0;
      chrisD[2][1][2][3] = 0.0;
      chrisD[2][1][3][0] = 0.0;
      chrisD[2][1][3][1] = 0.0;
      chrisD[2][1][3][2] = 0.0;
      chrisD[2][1][3][3] = 0.0;
      chrisD[2][2][0][0] = t14*p2*(2.0*p2-1.0)*t2;
      chrisD[2][2][0][1] = 0.0;
      chrisD[2][2][0][2] = 0.0;
      chrisD[2][2][0][3] = 0.0;
      chrisD[2][2][1][0] = 0.0;
      chrisD[2][2][1][1] = 0.0;
      chrisD[2][2][1][2] = 0.0;
      chrisD[2][2][1][3] = 0.0;
      chrisD[2][2][2][0] = 0.0;
      chrisD[2][2][2][1] = 0.0;
      chrisD[2][2][2][2] = 0.0;
      chrisD[2][2][2][3] = 0.0;
      chrisD[2][2][3][0] = 0.0;
      chrisD[2][2][3][1] = 0.0;
      chrisD[2][2][3][2] = 0.0;
      chrisD[2][2][3][3] = 0.0;
      chrisD[2][3][0][0] = 0.0;
      chrisD[2][3][0][1] = 0.0;
      chrisD[2][3][0][2] = 0.0;
      chrisD[2][3][0][3] = 0.0;
      chrisD[2][3][1][0] = 0.0;
      chrisD[2][3][1][1] = 0.0;
      chrisD[2][3][1][2] = 0.0;
      chrisD[2][3][1][3] = 0.0;
      chrisD[2][3][2][0] = 0.0;
      chrisD[2][3][2][1] = 0.0;
      chrisD[2][3][2][2] = 0.0;
      chrisD[2][3][2][3] = 0.0;
      chrisD[2][3][3][0] = 0.0;
      chrisD[2][3][3][1] = 0.0;
      chrisD[2][3][3][2] = 0.0;
      chrisD[2][3][3][3] = 0.0;
      chrisD[3][0][0][0] = 0.0;
      chrisD[3][0][0][1] = 0.0;
      chrisD[3][0][0][2] = 0.0;
      chrisD[3][0][0][3] = 0.0;
      chrisD[3][0][1][0] = 0.0;
      chrisD[3][0][1][1] = 0.0;
      chrisD[3][0][1][2] = 0.0;
      chrisD[3][0][1][3] = 0.0;
      chrisD[3][0][2][0] = 0.0;
      chrisD[3][0][2][1] = 0.0;
      chrisD[3][0][2][2] = 0.0;
      chrisD[3][0][2][3] = 0.0;
      chrisD[3][0][3][0] = -t5;
      chrisD[3][0][3][1] = 0.0;
      chrisD[3][0][3][2] = 0.0;
      chrisD[3][0][3][3] = 0.0;
      chrisD[3][1][0][0] = 0.0;
      chrisD[3][1][0][1] = 0.0;
      chrisD[3][1][0][2] = 0.0;
      chrisD[3][1][0][3] = 0.0;
      chrisD[3][1][1][0] = 0.0;
      chrisD[3][1][1][1] = 0.0;
      chrisD[3][1][1][2] = 0.0;
      chrisD[3][1][1][3] = 0.0;
      chrisD[3][1][2][0] = 0.0;
      chrisD[3][1][2][1] = 0.0;
      chrisD[3][1][2][2] = 0.0;
      chrisD[3][1][2][3] = 0.0;
      chrisD[3][1][3][0] = 0.0;
      chrisD[3][1][3][1] = 0.0;
      chrisD[3][1][3][2] = 0.0;
      chrisD[3][1][3][3] = 0.0;
      chrisD[3][2][0][0] = 0.0;
      chrisD[3][2][0][1] = 0.0;
      chrisD[3][2][0][2] = 0.0;
      chrisD[3][2][0][3] = 0.0;
      chrisD[3][2][1][0] = 0.0;
      chrisD[3][2][1][1] = 0.0;
      chrisD[3][2][1][2] = 0.0;
      chrisD[3][2][1][3] = 0.0;
      chrisD[3][2][2][0] = 0.0;
      chrisD[3][2][2][1] = 0.0;
      chrisD[3][2][2][2] = 0.0;
      chrisD[3][2][2][3] = 0.0;
      chrisD[3][2][3][0] = 0.0;
      chrisD[3][2][3][1] = 0.0;
      chrisD[3][2][3][2] = 0.0;
      chrisD[3][2][3][3] = 0.0;
      chrisD[3][3][0][0] = t21*p3*(2.0*p3-1.0)*t2;
      chrisD[3][3][0][1] = 0.0;
      chrisD[3][3][0][2] = 0.0;
      chrisD[3][3][0][3] = 0.0;
      chrisD[3][3][1][0] = 0.0;
      chrisD[3][3][1][1] = 0.0;
      chrisD[3][3][1][2] = 0.0;
      chrisD[3][3][1][3] = 0.0;
      chrisD[3][3][2][0] = 0.0;
      chrisD[3][3][2][1] = 0.0;
      chrisD[3][3][2][2] = 0.0;
      chrisD[3][3][2][3] = 0.0;
      chrisD[3][3][3][0] = 0.0;
      chrisD[3][3][3][1] = 0.0;
      chrisD[3][3][3][2] = 0.0;
      chrisD[3][3][3][3] = 0.0;
      
  return true;
}

/*! Transform local 4-direction to coordinate 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  ldir :  pointer to local direction array.
 *  \param  dir  :  pointer to calculated coordinate direction array.
 *  \param  type :  type of tetrad.
 */
void
MetricKasner :: localToCoord ( const double* pos, const double* ldir, double* dir, 
                               enum_nat_tetrad_type   )
{ 
  double t = pos[0];

  dir[0] = ldir[0];
  dir[1] = ldir[1]*pow(t,-p1);
  dir[2] = ldir[2]*pow(t,-p2);
  dir[3] = ldir[3]*pow(t,-p3);
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void 
MetricKasner :: coordToLocal ( const double* pos, const double* cdir, double* ldir, 
                               enum_nat_tetrad_type   )
{
  double t = pos[0];

  ldir[0] = cdir[0];
  ldir[1] = cdir[1]*pow(t,p1);
  ldir[2] = cdir[2]*pow(t,p2);
  ldir[3] = cdir[3]*pow(t,p3);
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool 
MetricKasner :: breakCondition ( const double* pos)
{
  bool br = false;

  double t = pos[0];
  if ((t<=0.0)) { br=true; }
  return br;
}


/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool 
MetricKasner :: setParam ( std::string pName, double val )
{
  Metric::setParam(pName,val);

  if (pName=="u")
  {
    mU = val;    
    calc_parameters();
  }
  return true;
}

/*! Generate report.
 */
bool    
MetricKasner :: report ( const vec4 , const vec4 , std::string &text )
{
  std::stringstream ss;
  ss << "Report for Kasner metric\n\tcoordinate : (t,x,y,z)\n";
  ss << "---------------------------------------------------------------\n";
  ss << "  physical units ..................... no\n";
  ss.precision(DEF_FIXED_REPORT_PRECISION);
  ss.setf(std::ios::fixed);
  ss << "  param .............................. u   = " << mU << std::endl;
  ss << "  param .............................. p_1 = " << p1 << std::endl;
  ss << "  param .............................. p_2 = " << p2 << std::endl;
  ss << "  param .............................. p_3 = " << p3 << std::endl;
  
  text = ss.str();
  return true;
}

/*! Calculate parameters p1,p2,p3 out of Khalatnikov-Lifshitz parameter u.
 */
void
MetricKasner :: calc_parameters ( )
{
   double edHN = 1.0/(1.0+mU+mU*mU);
   p1 = -mU*edHN;
   p2 = (1.0+mU)*edHN;
   p3 = mU*p2;
}

// ********************************* protected methods *****************************
/*!
 */
void 
MetricKasner :: setStandardValues( )
{
  mInitPos[0] = 1.0; 
  mInitPos[1] = 0.0;
  mInitPos[2] = 0.0;
  mInitPos[3] = 0.0;
  mInitDir[0] = 1.0;
  mInitDir[1] = 0.0;
  mInitDir[2] = 0.0;
   
  mCoordNames[0] = std::string("t");
  mCoordNames[1] = std::string("x");
  mCoordNames[2] = std::string("y");
  mCoordNames[3] = std::string("z"); 
}

} // end namespace m4d
