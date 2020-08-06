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
   m4dMetricMinkRotLattice.cpp
 
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

#include "m4dMetricMinkRotLattice.h"

namespace m4d
{
  
/*! Standard constructor for the rotating lattice metric.
 *
 * \param  omega : rotation parameter.
 */
MetricMinkRotLattice :: MetricMinkRotLattice( double omega )
{
  mMetricName = "RotLattice";
  setCoordType(enum_coordinate_cylinder);
  
  mOmega = omega;
  
  mPhysicalUnits = enum_physical_constants_geom;
  mSpeedOfLight = 1.0;
  mGravConstant = 1.0;
  
  mLocTeds.push_back(enum_nat_tetrad_static);
  mLocTeds.push_back(enum_nat_tetrad_lnrf);

  addParam("omega",mOmega);
  
  setStandardValues();
}

MetricMinkRotLattice :: ~MetricMinkRotLattice()
{
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricMinkRotLattice :: calculateMetric ( const double* pos )
{
  double c = mSpeedOfLight;
  double r = pos[1];
  
  double t1 = c*c;
  double t2 = r*r;
  double t3 = mOmega*mOmega;
  double t6 = mOmega*t2;
  
  g_compts[0][0] = -t1+t2*t3;
  g_compts[0][1] = 0.0;
  g_compts[0][2] = t6;
  g_compts[0][3] = 0.0;
  g_compts[1][0] = 0.0;
  g_compts[1][1] = 1.0;
  g_compts[1][2] = 0.0;
  g_compts[1][3] = 0.0;
  g_compts[2][0] = t6;
  g_compts[2][1] = 0.0;
  g_compts[2][2] = t2;
  g_compts[2][3] = 0.0;
  g_compts[3][0] = 0.0;
  g_compts[3][1] = 0.0;
  g_compts[3][2] = 0.0;
  g_compts[3][3] = 1.0;

  return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricMinkRotLattice :: calculateChristoffels  ( const double* pos )
{
  double r = pos[1];
  
  double t1 = mOmega*mOmega;
  double t3 = 1/r;
  double t4 = mOmega*t3;
  double t5 = mOmega*r;
  
  christoffel[0][0][0] = 0.0;
  christoffel[0][0][1] = -r*t1;
  christoffel[0][0][2] = 0.0;
  christoffel[0][0][3] = 0.0;
  christoffel[0][1][0] = 0.0;
  christoffel[0][1][1] = 0.0;
  christoffel[0][1][2] = t4;
  christoffel[0][1][3] = 0.0;
  christoffel[0][2][0] = 0.0;
  christoffel[0][2][1] = -t5;
  christoffel[0][2][2] = 0.0;
  christoffel[0][2][3] = 0.0;
  christoffel[0][3][0] = 0.0;
  christoffel[0][3][1] = 0.0;
  christoffel[0][3][2] = 0.0;
  christoffel[0][3][3] = 0.0;
  christoffel[1][0][0] = 0.0;
  christoffel[1][0][1] = 0.0;
  christoffel[1][0][2] = t4;
  christoffel[1][0][3] = 0.0;
  christoffel[1][1][0] = 0.0;
  christoffel[1][1][1] = 0.0;
  christoffel[1][1][2] = 0.0;
  christoffel[1][1][3] = 0.0;
  christoffel[1][2][0] = 0.0;
  christoffel[1][2][1] = 0.0;
  christoffel[1][2][2] = t3;
  christoffel[1][2][3] = 0.0;
  christoffel[1][3][0] = 0.0;
  christoffel[1][3][1] = 0.0;
  christoffel[1][3][2] = 0.0;
  christoffel[1][3][3] = 0.0;
  christoffel[2][0][0] = 0.0;
  christoffel[2][0][1] = -t5;
  christoffel[2][0][2] = 0.0;
  christoffel[2][0][3] = 0.0;
  christoffel[2][1][0] = 0.0;
  christoffel[2][1][1] = 0.0;
  christoffel[2][1][2] = t3;
  christoffel[2][1][3] = 0.0;
  christoffel[2][2][0] = 0.0;
  christoffel[2][2][1] = -r;
  christoffel[2][2][2] = 0.0;
  christoffel[2][2][3] = 0.0;
  christoffel[2][3][0] = 0.0;
  christoffel[2][3][1] = 0.0;
  christoffel[2][3][2] = 0.0;
  christoffel[2][3][3] = 0.0;
  christoffel[3][0][0] = 0.0;
  christoffel[3][0][1] = 0.0;
  christoffel[3][0][2] = 0.0;
  christoffel[3][0][3] = 0.0;
  christoffel[3][1][0] = 0.0;
  christoffel[3][1][1] = 0.0;
  christoffel[3][1][2] = 0.0;
  christoffel[3][1][3] = 0.0;
  christoffel[3][2][0] = 0.0;
  christoffel[3][2][1] = 0.0;
  christoffel[3][2][2] = 0.0;
  christoffel[3][2][3] = 0.0;
  christoffel[3][3][0] = 0.0;
  christoffel[3][3][1] = 0.0;
  christoffel[3][3][2] = 0.0;
  christoffel[3][3][3] = 0.0;
  
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
MetricMinkRotLattice :: localToCoord ( const double* pos, const double* ldir, double* dir,
                                   enum_nat_tetrad_type  type  )
{
  double r = pos[1];
  double wdc = mOmega/mSpeedOfLight;
  
  if (type == enum_nat_tetrad_default || type == enum_nat_tetrad_static)
  {
    double sq = sqrt(1.0-r*r*wdc*wdc);
    
    dir[0] = ldir[0]/sq/mSpeedOfLight + ldir[2]*wdc*r/sq/mSpeedOfLight;
    dir[1] = ldir[1];
    dir[2] = ldir[2]*sq/r;
    dir[3] = ldir[3];
  }
  else if (type == enum_nat_tetrad_lnrf)
  {  
    dir[0] = ldir[0]/mSpeedOfLight;
    dir[1] = ldir[1];
    dir[2] = ldir[2]/r-wdc*ldir[0];
    dir[3] = ldir[3];
  }
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void
MetricMinkRotLattice :: coordToLocal ( const double* pos, const double* cdir, double* ldir,
                                   enum_nat_tetrad_type  type  )
{
  double r = pos[1];
  double wdc = mOmega/mSpeedOfLight;
  
  if (type == enum_nat_tetrad_default || type == enum_nat_tetrad_static)
  {
    double sq = sqrt(1.0-r*r*wdc*wdc);
    
    ldir[0] = cdir[0]*mSpeedOfLight*sq - cdir[2]*wdc*r*r/sq;
    ldir[1] = cdir[1];
    ldir[2] = r/sq*cdir[2];
    ldir[3] = cdir[3];
  }
  else if (type == enum_nat_tetrad_lnrf)
  {
    ldir[0] = cdir[0]*mSpeedOfLight;
    ldir[1] = cdir[1];
    ldir[2] = r*(cdir[2] + mOmega*cdir[0]);
    ldir[3] = cdir[3];
  }
}

/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position > c/omega.
 *  \return false : position is valid.
 */
bool 
MetricMinkRotLattice :: breakCondition ( const double* pos)
{  
  if (pos[1] >= fabs(mSpeedOfLight/mOmega))
    return true;

  return false;
}



/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'omega' parameter.
 */
bool 
MetricMinkRotLattice :: setParam ( std::string pName, double val )
{
  if (Metric::setParam(pName,val))
  {
    mOmega = val;
  }
  return true;
}

/*! Tests whether the constraint equation is fulfilled.
 *
 *  The constraint equation for lightlike and timelike geodesics reads:  
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  However, take care of the limited double precision.
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double 
MetricMinkRotLattice :: testConstraint ( const double y[], const double kappa )
{
  double r = y[1];
  double cm = 1.0/mSpeedOfLight;
  
  // Scale the directions with the speed of light before doubling them !!
  double dt  = y[4];
  double dr  = y[5]*cm;
  double dph = y[6]*cm;
  double dz  = y[7]*cm;
  
  double wdc = mOmega*cm;

  double sum = -kappa;  
  sum += -(1.0-wdc*wdc*r*r)*dt*dt + 2.0*r*r*cm*mOmega*dt*dph + r*r*dph*dph +  dr*dr + dz*dz;
  
  return sum;
}

/*! Generate report.
 */
bool    
MetricMinkRotLattice :: report ( const vec4 , const vec4 , std::string &text )
{
  std::stringstream ss;
  ss << "Report for Minkowski as rotating lattice metric\n\tcoordinate : (t,r,phi,z)\n";
  ss << "--------------------------------------------------------\n";
  ss << "  physical units ................................. no\n";
  ss << "  points outside  r=c/omega are not allowed ...... " << mSpeedOfLight/mOmega << std::endl;
  text = ss.str();
  return true;
}

// ********************************* protected methods *****************************
/*!
 */
void
MetricMinkRotLattice :: setStandardValues ( )
{

  mInitPos[0] = 0.0; 
  mInitPos[1] = 1.0;
  mInitPos[2] = 0.0;
  mInitPos[3] = 0.0;
  mInitDir[0] = 1.0;
  mInitDir[1] = 0.0;
  mInitDir[2] = 0.0;
   
  mCoordNames[0] = std::string("t");
  mCoordNames[1] = std::string("r");
  mCoordNames[2] = std::string("phi");
  mCoordNames[3] = std::string("z"); 
}

} // end namespace m4d
