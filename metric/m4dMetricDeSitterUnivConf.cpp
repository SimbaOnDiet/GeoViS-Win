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
   m4dMetricDeSitterUnivConf.cpp
 
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

#include "m4dMetricDeSitterUnivConf.h"

namespace m4d
{

#define eps 1.0e-6


/*! Standard constructor for the Kottler metric.
 *
 * \param  h : Hubble parameter.
 * \param  p : artificial parameter to switch the conformal factor.
 */
MetricDeSitterUnivConf :: MetricDeSitterUnivConf ( double h, double p )
{
  mMetricName  = "DeSitterUnivConformal";
  setCoordType(enum_coordinate_cartesian);
    
  mPhysicalUnits = enum_physical_constants_geom;
  mSpeedOfLight = 1.0;
  mGravConstant = 1.0;
  
  addParam("hubble",h);
  mHubble = h;

  addParam("p",p);
  mP = p;
  
  mLocTeds.push_back(enum_nat_tetrad_static);

  mDrawTypes.push_back(enum_draw_twoplusone);

  setStandardValues();
}

MetricDeSitterUnivConf :: ~MetricDeSitterUnivConf()
{
}


// *********************************** public methods ******************************
/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricDeSitterUnivConf :: calculateMetric  ( const double* pos )
{
  double T = pos[0];
	
  double H = mHubble;
  double c = mSpeedOfLight;	
	  
  double t1 = c*c;
  double t2 = H*H;
  double t5 = T*T;
  double t7 = t1/t2/t5;
	
  if (mP<0.0)
  {
	t7 = 1.0;  
  }
  
      g_compts[0][0] = -t7;
      g_compts[0][1] = 0.0;
      g_compts[0][2] = 0.0;
      g_compts[0][3] = 0.0;
      g_compts[1][0] = 0.0;
      g_compts[1][1] = t7;
      g_compts[1][2] = 0.0;
      g_compts[1][3] = 0.0;
      g_compts[2][0] = 0.0;
      g_compts[2][1] = 0.0;
      g_compts[2][2] = t7;
      g_compts[2][3] = 0.0;
      g_compts[3][0] = 0.0;
      g_compts[3][1] = 0.0;
      g_compts[3][2] = 0.0;
      g_compts[3][3] = t7;	  

  return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricDeSitterUnivConf :: calculateChristoffels  ( const double* pos )
{
  double T = pos[0];
  double t1 = 1/T;
  
  if (mP<0.0)
  {
	  t1 = 0.0;
  }
  
      christoffel[0][0][0] = -t1;
      christoffel[0][0][1] = 0.0;
      christoffel[0][0][2] = 0.0;
      christoffel[0][0][3] = 0.0;
      christoffel[0][1][0] = 0.0;
      christoffel[0][1][1] = -t1;
      christoffel[0][1][2] = 0.0;
      christoffel[0][1][3] = 0.0;
      christoffel[0][2][0] = 0.0;
      christoffel[0][2][1] = 0.0;
      christoffel[0][2][2] = -t1;
      christoffel[0][2][3] = 0.0;
      christoffel[0][3][0] = 0.0;
      christoffel[0][3][1] = 0.0;
      christoffel[0][3][2] = 0.0;
      christoffel[0][3][3] = -t1;
      christoffel[1][0][0] = 0.0;
      christoffel[1][0][1] = -t1;
      christoffel[1][0][2] = 0.0;
      christoffel[1][0][3] = 0.0;
      christoffel[1][1][0] = -t1;
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
      christoffel[2][0][2] = -t1;
      christoffel[2][0][3] = 0.0;
      christoffel[2][1][0] = 0.0;
      christoffel[2][1][1] = 0.0;
      christoffel[2][1][2] = 0.0;
      christoffel[2][1][3] = 0.0;
      christoffel[2][2][0] = -t1;
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
      christoffel[3][0][3] = -t1;
      christoffel[3][1][0] = 0.0;
      christoffel[3][1][1] = 0.0;
      christoffel[3][1][2] = 0.0;
      christoffel[3][1][3] = 0.0;
      christoffel[3][2][0] = 0.0;
      christoffel[3][2][1] = 0.0;
      christoffel[3][2][2] = 0.0;
      christoffel[3][2][3] = 0.0;
      christoffel[3][3][0] = -t1;
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
MetricDeSitterUnivConf :: calculateChrisD ( const double* pos )
{
  double T = pos[0];
  
  double t1 = T*T;
  double t2 = 1/t1;

  if (mP<0.0)
  {
	  t1 = t2 = 0.0;
  }
      chrisD[0][0][0][0] = t2;
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
      chrisD[0][1][1][0] = t2;
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
      chrisD[0][2][2][0] = t2;
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
      chrisD[0][3][3][0] = t2;
      chrisD[0][3][3][1] = 0.0;
      chrisD[0][3][3][2] = 0.0;
      chrisD[0][3][3][3] = 0.0;
      chrisD[1][0][0][0] = 0.0;
      chrisD[1][0][0][1] = 0.0;
      chrisD[1][0][0][2] = 0.0;
      chrisD[1][0][0][3] = 0.0;
      chrisD[1][0][1][0] = t2;
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
      chrisD[1][1][0][0] = t2;
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
      chrisD[2][0][2][0] = t2;
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
      chrisD[2][2][0][0] = t2;
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
      chrisD[3][0][3][0] = t2;
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
      chrisD[3][3][0][0] = t2;
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
MetricDeSitterUnivConf :: localToCoord ( const double* pos, const double* ldir, double* dir, 
                                         enum_nat_tetrad_type   )
{ 
  double T = pos[0];  
  double f = -mHubble*T/mSpeedOfLight;
  
  if (mP<0.0)
    f = 1.0;
  
  dir[0] = ldir[0]*f;
  dir[1] = ldir[1]*f;
  dir[2] = ldir[2]*f;
  dir[3] = ldir[3]*f;  
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void 
MetricDeSitterUnivConf :: coordToLocal ( const double* pos, const double* cdir, double* ldir, 
                                         enum_nat_tetrad_type   )
{
  double T = pos[0];  
  double edf = -mSpeedOfLight/(mHubble*T);
  
  if (mP<0.0)
    edf = 1.0;
  
  ldir[0] = cdir[0]*edf;
  ldir[1] = cdir[1]*edf;
  ldir[2] = cdir[2]*edf;
  ldir[3] = cdir[3]*edf;
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool 
MetricDeSitterUnivConf :: breakCondition ( const double* pos)
{
  bool br = false;

  double T = pos[0];
  

  if ((T+1e-12)>=0.0 && !(mP<0.0)) br=true;
  return br;
}

/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:  
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double 
MetricDeSitterUnivConf :: testConstraint ( const double y[], const double kappa )
{
  double T = y[0];
  double edf = mSpeedOfLight/(mHubble*T);
  double cm  = 1.0/mSpeedOfLight;

  if (mP<0.0)
    edf = 1.0;
    
  // Scale the directions with the speed of light before doubling them !!
  double dT = y[4];
  double dx = y[5]*cm;
  double dy = y[6]*cm;
  double dz = y[7]*cm;

  double sum = -kappa;
  sum += edf*edf*(-dT*dT + dx*dx + dy*dy + dz*dz);
  return sum;
}


/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool
MetricDeSitterUnivConf :: transToTwoPlusOne ( vec4 p, vec4 &cp )
{
  cp = vec4( p[0], p[1], p[2], p[0] );
  return true;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' or 'lambda' parameter.
 */
bool 
MetricDeSitterUnivConf :: setParam ( std::string pName, double val )
{
  Metric::setParam(pName,val);

  if (pName=="hubble")
    mHubble = val;
  else if (pName=="p")
    mP = val;
  return true;
}


/*! Generate report.
 */
bool    
MetricDeSitterUnivConf :: report ( const vec4 , const vec4 , std::string &text )
{
  std::stringstream ss;
  ss << "Report for conformal deSitter universe metric\n\tcoordinate : (T,x,y,z)\n";
  ss << "---------------------------------------------------------\n";
  ss << "  physical units ................................. no\n";
  
  text = ss.str();
  return true;
}


// ********************************* protected methods *****************************
/*!
 */
void 
MetricDeSitterUnivConf :: setStandardValues( )
{
  mInitPos[0] = -1.0; 
  mInitPos[1] = 0.0;
  mInitPos[2] = 0.0;
  mInitPos[3] = 0.0;
  mInitDir[0] = 1.0;
  mInitDir[1] = 0.0;
  mInitDir[2] = 0.0;
   
  mCoordNames[0] = std::string("T");
  mCoordNames[1] = std::string("x");
  mCoordNames[2] = std::string("y");
  mCoordNames[3] = std::string("z"); 
}

} // end namespace m4d
