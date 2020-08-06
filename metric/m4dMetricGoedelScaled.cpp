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
    m4dMetricGoedelScaled.cpp 
 
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

#include "m4dMetricGoedelScaled.h"

namespace m4d
{
//   
// ---------------------------------------------------
//    constructur/destructor
// ---------------------------------------------------
MetricGoedelScaled :: MetricGoedelScaled ( double rG, double zeta )
{
  mMetricName  = "Goedel scaled";
  setCoordType(enum_coordinate_cylinder);
    
  mPhysicalUnits = enum_physical_constants_geom;
  mSpeedOfLight = 1.0;
  mGravConstant = 1.0;
  
  addParam("rG");
  setParam("rG",rG);
  addParam("zeta");
  setParam("zeta",zeta);
  
  mRG = rG;
  mZeta = zeta;

  mDrawTypes.push_back(enum_draw_twoplusone);
  
  setStandardValues();
}

MetricGoedelScaled :: ~MetricGoedelScaled()
{
 
}


// *********************************** public methods ******************************
// ---------------------------------------------------
//    public :: calculateMetric
// ---------------------------------------------------
bool MetricGoedelScaled :: calculateMetric  ( const double* pos )
{
  double r_G = mRG;
  double R   = pos[1];
  double c   = mSpeedOfLight;

  double t1 = r_G*r_G;
  double t2 = c*c;
  double t4 = sqrt(2.0);
  double t6 = R*R;
  double t8 = t1*t4*c*t6;
  double t13 = t6*t6;
  
  g_compts[0][0] = -t1*t2;
  g_compts[0][1] = 0.0;
  g_compts[0][2] = -t8;
  g_compts[0][3] = 0.0;
  g_compts[1][0] = 0.0;
  g_compts[1][1] = t1/(1.0+t6);
  g_compts[1][2] = 0.0;
  g_compts[1][3] = 0.0;
  g_compts[2][0] = -t8;
  g_compts[2][1] = 0.0;
  g_compts[2][2] = t1*t6-t1*t13;
  g_compts[2][3] = 0.0;
  g_compts[3][0] = 0.0;
  g_compts[3][1] = 0.0;
  g_compts[3][2] = 0.0;
  g_compts[3][3] = t1;

  return true; 
}

// ---------------------------------------------------
//    public :: calculateChristoffels
// ---------------------------------------------------
bool MetricGoedelScaled :: calculateChristoffels  ( const double* pos )
{
  double R = pos[1];
  double c = mSpeedOfLight;

  double t1 = R*R;
  double t2 = t1+1.0;
  double t3 = 1.0/t2;
  double t4 = t3*R;
  double t5 = 2.0*t4;
  double t7 = 1.0/R*t3;
  double t8 = sqrt(2.0);
  double t10 = t7*t8*c;
  double t13 = t2*t8*c*R;
  double t18 = t8*t1*R/c*t3;

  christoffel[0][0][0] = 0.0;
  christoffel[0][0][1] = 0.0;
  christoffel[0][0][2] = 0.0;
  christoffel[0][0][3] = 0.0;
  christoffel[0][1][0] = t5;
  christoffel[0][1][1] = 0.0;
  christoffel[0][1][2] = -t10;
  christoffel[0][1][3] = 0.0;
  christoffel[0][2][0] = 0.0;
  christoffel[0][2][1] = t13;
  christoffel[0][2][2] = 0.0;
  christoffel[0][2][3] = 0.0;
  christoffel[0][3][0] = 0.0;
  christoffel[0][3][1] = 0.0;
  christoffel[0][3][2] = 0.0;
  christoffel[0][3][3] = 0.0;
  christoffel[1][0][0] = t5;
  christoffel[1][0][1] = 0.0;
  christoffel[1][0][2] = -t10;
  christoffel[1][0][3] = 0.0;
  christoffel[1][1][0] = 0.0;
  christoffel[1][1][1] = -t4;
  christoffel[1][1][2] = 0.0;
  christoffel[1][1][3] = 0.0;
  christoffel[1][2][0] = t18;
  christoffel[1][2][1] = 0.0;
  christoffel[1][2][2] = t7;
  christoffel[1][2][3] = 0.0;
  christoffel[1][3][0] = 0.0;
  christoffel[1][3][1] = 0.0;
  christoffel[1][3][2] = 0.0;
  christoffel[1][3][3] = 0.0;
  christoffel[2][0][0] = 0.0;
  christoffel[2][0][1] = t13;
  christoffel[2][0][2] = 0.0;
  christoffel[2][0][3] = 0.0;
  christoffel[2][1][0] = t18;
  christoffel[2][1][1] = 0.0;
  christoffel[2][1][2] = t7;
  christoffel[2][1][3] = 0.0;
  christoffel[2][2][0] = 0.0;
  christoffel[2][2][1] = t2*R*(-1.0+2.0*t1);
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


// ---------------------------------------------------
//    public :: localToCoord
// ---------------------------------------------------
void MetricGoedelScaled :: localToCoord   ( const double* pos, const double* ldir, double* dir,
                                               enum_nat_tetrad_type   )
{
  double R  = pos[1];
  double R2 = R*R;
  double c  = mSpeedOfLight;
  double c2 = c*c;

  double A = R2*(-sqrt(2.0)*c+(1.0-R2)*mZeta);
  double B = c2 + sqrt(2.0)*R2*c*mZeta;
  double Gamma = 1.0/sqrt(c2 + 2.0*sqrt(2.0)*R2*c*mZeta - R2*(1.0-R2)*mZeta*mZeta );
  double Delta = 1.0/(R*c*sqrt(1.0+R2));

  dir[0] = Gamma/mRG*ldir[0]       + Delta*Gamma*A/mRG*ldir[2];
  dir[1] = 1.0/mRG*sqrt(1.0+R2)*ldir[1];
  dir[2] = Gamma*mZeta/mRG*ldir[0] + Delta*Gamma*B/mRG*ldir[2];
  dir[3] = 1.0/mRG*ldir[3];

}

// ---------------------------------------------------
//    public :: coordToLocal
// ---------------------------------------------------
void MetricGoedelScaled :: coordToLocal   ( const double* pos, const double* cdir, double* ldir,
                                               enum_nat_tetrad_type   )
{
  double R  = pos[1];
  double R2 = R*R;
  double c  = mSpeedOfLight;
  double c2 = c*c;

  double A = R2*(-sqrt(2.0)*c+(1.0-R2)*mZeta);
  double B = c2 + sqrt(2.0)*R2*c*mZeta;
  double Gamma = 1.0/sqrt(c2 + 2.0*sqrt(2.0)*R2*c*mZeta - R2*(1.0-R2)*mZeta*mZeta );
  double Delta = 1.0/(R*c*sqrt(1.0+R2));

  double nenner = B-mZeta*A;

  ldir[0] = mRG/Gamma         * (B*cdir[0] - A*cdir[2]   )/nenner;
  ldir[1] = mRG/sqrt(1.0+R2)*cdir[1];
  ldir[2] = mRG/(Gamma*Delta) * (cdir[2] - mZeta*cdir[0] )/nenner;
  ldir[3] = mRG*cdir[3];
}


// ---------------------------------------------------
//    public :: breakCondition
// ---------------------------------------------------
bool MetricGoedelScaled :: breakCondition ( const double* )
{
  return false;
}

// ---------------------------------------------------
//    public :: setParam
// ---------------------------------------------------
bool MetricGoedelScaled :: setParam ( std::string pName, double val )
{
  Metric::setParam(pName,val);
  if (pName == "rG")
    mRG = val;
  else if (pName == "zeta")
    mZeta = val;

  return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool
MetricGoedelScaled :: transToTwoPlusOne  ( vec4 p, vec4 &cp )
{
  vec4 tp;
  TransCoordinates::toCartesianCoord( mCoordType, p, tp );
  cp = vec4( tp[0], tp[1], tp[2], tp[0] );
  return true;
}



/*! Generate report.
 */
bool
MetricGoedelScaled :: report ( const vec4 pos, const vec4 cdir, std::string &text )
{
  std::stringstream ss;
  ss << "Report for Goedel metric\n\tcoordinate : (T,R,Phi,Z)\n";
  ss << "\tScaled coordinates, geodesical shape independent of rG\n";
  ss << "---------------------------------------------------------------\n";
  ss << "  physical units ........... yes\n";
  ss.precision(DEF_FIXED_REPORT_PRECISION);
  ss.setf(std::ios::fixed);
  ss.precision(DEF_FIXED_REPORT_PRECISION);
  ss.setf(std::ios::fixed);
  vec4 locStartDir;
  coordToLocal(pos.data(),cdir.data(),locStartDir.data());
  double k0  = -locStartDir[0]/mRG;
  double k2  = (pos[1]*sqrt(1.0+pos[1]*pos[1])*locStartDir[2] - sqrt(2.0)*pos[1]*pos[1]*locStartDir[0])/mRG;
  double k3  = locStartDir[3]/mRG;
  ss << "  Goedel radius ............ rG = 2a = " << mRG << std::endl;
  ss << "  constant of motion ....... k_0 = " << k0 << std::endl;
  ss << "  constant of motion ....... k_2 = " << k2 << std::endl;
  ss << "  constant of motion ....... k_3 = " << k3;

  text = ss.str();
  return true;
}


// ********************************* protected methods *****************************
// ---------------------------------------------------
//    protected :: setViewerVal
// ---------------------------------------------------
void 
MetricGoedelScaled :: setStandardValues ( )
{
  mInitPos[0] = 0.0; 
  mInitPos[1] = 0.1;
  mInitPos[2] = 0.0;
  mInitPos[3] = 0.0;
  mInitDir[0] = 1.0;
  mInitDir[1] = 0.0;
  mInitDir[2] = 0.0;
   
  mCoordNames[0] = std::string("T");
  mCoordNames[1] = std::string("R");
  mCoordNames[2] = std::string("phi");
  mCoordNames[3] = std::string("Z");
}

} // end namespace m4d
