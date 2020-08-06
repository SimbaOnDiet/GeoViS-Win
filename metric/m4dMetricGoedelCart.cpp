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
    m4dMetricGoedelCart.cpp
 
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

#include "m4dMetricGoedelCart.h"

namespace m4d
{
  
// ---------------------------------------------------
//    constructur/destructor
// ---------------------------------------------------
MetricGoedelCart :: MetricGoedelCart ( double a, double zeta )
{
  mMetricName  = "Goedel cartesian";
  setCoordType(enum_coordinate_cartesian);
    
  mPhysicalUnits = enum_physical_constants_geom;
  mSpeedOfLight = 1.0;
  mGravConstant = 1.0;
  
  addParam("a");
  setParam("a",a);
  addParam("zeta");
  setParam("zeta",zeta);
  
  mDrawTypes.push_back(enum_draw_twoplusone);

  setStandardValues();
}

MetricGoedelCart :: ~MetricGoedelCart()
{
  
}


// *********************************** public methods ******************************
// ---------------------------------------------------
//    public :: calculateMetric
// ---------------------------------------------------
bool MetricGoedelCart
  :: calculateMetric  ( const double* pos )
{
  double x = pos[1];
  double y = pos[2];

  double x2 = x*x;
  double y2 = y*y;
  double r = sqrt(x*x+y*y);
  double r2 = r*r;
  
  double r2a  = r/(2.0*mA);
  double fac1 = 1.0/(1.0+r2a*r2a);
  double fac2 = 1.0-r2a*r2a;   
    
  g_compts[0][0] = -mSpeedOfLight*mSpeedOfLight;
  g_compts[0][1] =  mSpeedOfLight/(sqrt(2.0)*mA) * y;
  g_compts[0][2] = -mSpeedOfLight/(sqrt(2.0)*mA) * x;
  g_compts[0][3] = 0.0;
  g_compts[1][0] =  mSpeedOfLight/(sqrt(2.0)*mA) * y;
  g_compts[1][1] = (fac1*x2 + fac2*y2)/r2;
  g_compts[1][2] = (fac1-fac2)*x*y/r2;
  g_compts[1][3] = 0.0;
  g_compts[2][0] = -mSpeedOfLight/(sqrt(2.0)*mA) * x;
  g_compts[2][1] = (fac1-fac2)*x*y/r2;
  g_compts[2][2] = (fac1*y2 + fac2*x2)/r2;
  g_compts[2][3] = 0.0;
  g_compts[3][0] = 0.0;
  g_compts[3][1] = 0.0;
  g_compts[3][2] = 0.0;
  g_compts[3][3] = 1.0;      
  
  return true;
}

// ---------------------------------------------------
//    public :: calculateChristoffels
// ---------------------------------------------------
bool MetricGoedelCart
  :: calculateChristoffels  ( const double* pos )
{
  // no coordinate singularities in cartesian coordinates  
  double a = mA;
  double c = mSpeedOfLight;
  double x = pos[1];
  double y = pos[2];

  double t1 = a*a;
  double t3 = x*x;
  double t4 = y*y;
  double t6 = 1/(4.0*t1+t3+t4);
  double t8 = 2.0*x*t6;
  double t14 = 1/t1/a;
  double t16 = sqrt(2.0);
  double t20 = y*x*(8.0*t1+t3+t4)*t14*t6*c*t16/8.0;
  double t21 = t3*t4;
  double t22 = t4*t4;
  double t23 = t1*t1;
  double t24 = 16.0*t23;
  double t25 = t1*t4;
  double t30 = t6*c*t16;
  double t32 = (t21+t22+t24+8.0*t25)*t14*t30/8.0;
  double t34 = 2.0*y*t6;
  double t35 = t3*t3;
  double t36 = t1*t3;
  double t41 = (t35+t21+t24+8.0*t36)*t14*t30/8.0;
  double t47 = t6/a/c;
  double t48 = t16*y*x*t47;
  double t49 = 8.0*t23;
  double t50 = 6.0*t25;
  double t54 = 1/t23*t6;
  double t57 = t49+t50+t21+t22;
  double t64 = (t3-t4)*t16*t47/2.0;
  double t65 = 6.0*t36;
  double t66 = t49+t65+t35+t21;
  double t69 = t66*y*t54/8.0;
  double t72 = t57*x*t54/8.0;
  
  christoffel[0][0][0] = 0.0;
  christoffel[0][0][1] = 0.0;
  christoffel[0][0][2] = 0.0;
  christoffel[0][0][3] = 0.0;
  christoffel[0][1][0] = t8;
  christoffel[0][1][1] = -t20;
  christoffel[0][1][2] = -t32;
  christoffel[0][1][3] = 0.0;
  christoffel[0][2][0] = t34;
  christoffel[0][2][1] = t41;
  christoffel[0][2][2] = t20;
  christoffel[0][2][3] = 0.0;
  christoffel[0][3][0] = 0.0;
  christoffel[0][3][1] = 0.0;
  christoffel[0][3][2] = 0.0;
  christoffel[0][3][3] = 0.0;
  christoffel[1][0][0] = t8;
  christoffel[1][0][1] = -t20;
  christoffel[1][0][2] = -t32;
  christoffel[1][0][3] = 0.0;
  christoffel[1][1][0] = -t48;
  christoffel[1][1][1] = -(t49-t50-t21-t22)*x*t54/8.0;
  christoffel[1][1][2] = t57*y*t54/8.0;
  christoffel[1][1][3] = 0.0;
  christoffel[1][2][0] = t64;
  christoffel[1][2][1] = -t69;
  christoffel[1][2][2] = -t72;
  christoffel[1][2][3] = 0.0;
  christoffel[1][3][0] = 0.0;
  christoffel[1][3][1] = 0.0;
  christoffel[1][3][2] = 0.0;
  christoffel[1][3][3] = 0.0;
  christoffel[2][0][0] = t34;
  christoffel[2][0][1] = t41;
  christoffel[2][0][2] = t20;
  christoffel[2][0][3] = 0.0;
  christoffel[2][1][0] = t64;
  christoffel[2][1][1] = -t69;
  christoffel[2][1][2] = -t72;
  christoffel[2][1][3] = 0.0;
  christoffel[2][2][0] = t48;
  christoffel[2][2][1] = t66*x*t54/8.0;
  christoffel[2][2][2] = -(t49-t65-t35-t21)*y*t54/8.0;
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
void MetricGoedelCart
  :: localToCoord   ( const double* pos, const double* ldir, double* dir, 
                                               enum_nat_tetrad_type   )
{ 
  double x = pos[1];
  double y = pos[2];
  
 

  double r = sqrt(x*x+y*y);
  //double r2 = r*r;
  double r2a  = r/(2.0*mA);

  if ( r > m4dGoedelCartEps )
  {
    //full metric
    calcVars(pos);
    
    //local to coord to cylindrical coordinates + transformation to cartesian coordinate (direction)
    
    dir[0] = ldir[0]*Gamma      + ldir[2]*Gamma*Delta*F1;
    dir[1] = ldir[1]*sqrt(1.0+r2a*r2a)*x/r + (ldir[0]*Gamma*mZeta + ldir[2]*Gamma*Delta*F2)*(-y);
    dir[2] = ldir[1]*sqrt(1.0+r2a*r2a)*y/r + (ldir[0]*Gamma*mZeta + ldir[2]*Gamma*Delta*F2)*x;
    dir[3] = ldir[3]; 
    
//     printf("Gamma Delta F1 F2: %f %f %f %f\n",Gamma,Delta,F1,F2); 
  }
  else
  {
    //metric converges to minkowski metric   
    dir[0] = ldir[0];
    dir[1] = ldir[1];
    dir[2] = ldir[2];
    dir[3] = ldir[3]; 
  }  
 
}
// ---------------------------------------------------
//    public :: coordToLocal
// ---------------------------------------------------
void MetricGoedelCart
  :: coordToLocal   ( const double* , const double* , double* ,
                                               enum_nat_tetrad_type   )
{
  printf("MetricGoedelCart::coordToLocal missing\n");
}

// ---------------------------------------------------
//    public :: breakCondition
// ---------------------------------------------------
bool MetricGoedelCart
  :: breakCondition ( const double* )
{
  return false;
}

// ---------------------------------------------------
//    public :: setParam
// ---------------------------------------------------
bool MetricGoedelCart
  :: setParam ( std::string pName, double val )
{   
  Metric::setParam(pName,val);        
  if (pName == "a")
    mA = val;
  else if (pName == "zeta")
    mZeta = val;      
           
  //printf("a=%f, zeta=%f\n",mA,mZeta);
  
  return true;
}

/*! Transform point p to 2+1 coordinates.
 *
 *  \param  p  : point in proper metric coordinates.
 *  \param  cp : reference to transformed point.
 *  \return true : success.
 */
bool    
MetricGoedelCart :: transToTwoPlusOne  ( vec4 p, vec4 &cp )
{
  cp = vec4( p[0], p[1], p[2], p[0] );
  return true;
}

/*! Generate report.
 */
bool
MetricGoedelCart :: report ( const vec4 pos, const vec4 cdir, std::string &text )
{
  std::stringstream ss;
  ss << "Report for Goedel metric\n\tcoordinate : (t,x,y,z)\n";
  ss << "\tIncreasing parameter a -> straight geodesics\n";
  ss << "\talso valid at x=y=0 (no coordinate singularity\n";
  ss << "---------------------------------------------------------------\n";
  ss << "  physical units ........... yes\n";
  ss.precision(DEF_FIXED_REPORT_PRECISION);
  ss.setf(std::ios::fixed);
  vec4 locStartDir;
  coordToLocal(pos.data(),cdir.data(),locStartDir.data());
  double k0  = -locStartDir[0];
  double k2  = (pos[1]*sqrt(1.0+pos[1]*pos[1]/4.0/mA/mA)*locStartDir[2] - sqrt(2.0)*pos[1]*pos[1]*locStartDir[0])/2.0/mA;
  double k3  = locStartDir[3];
  ss << "  Goedel radius ............ rG = 2a = " << 2.0*mA << std::endl;
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
MetricGoedelCart :: setStandardValues ( )
{
  mInitPos[0] = 0.0; 
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

// ---------------------------------------------------
//    protected :: calcVars
// ---------------------------------------------------
void 
MetricGoedelCart :: calcVars ( const double* pos)
{
  double x = pos[1];
  double y = pos[2];

  double r = sqrt(x*x+y*y);
  double r2 = r*r;
  double r2a  = r/(2.0*mA);
  
  if ( r > m4dGoedelCartEps )
  {
    //full metric
    Gamma = sqrt(
                        + mSpeedOfLight*mSpeedOfLight 
                        + mZeta*r2*mSpeedOfLight*sqrt(2.0)/mA
                        - mZeta*mZeta*r2*(1.0-r2a*r2a)
                       );
    Gamma = 1.0/Gamma;
    
    Delta = r*mSpeedOfLight*sqrt(1.0+r2a*r2a);
    Delta = 1.0/Delta;
    
    F1 = ( -r2*mSpeedOfLight/(sqrt(2.0)*mA) + mZeta*r2*(1.0-r2a*r2a) );
    F2 = ( mSpeedOfLight*mSpeedOfLight + mZeta*r2*mSpeedOfLight/(sqrt(2.0)*mA) );
  }  
  else
  {
    printf("error in MetricGoedelCart :: setStandardValues -> r too small\n");
  }
  
  
}


} // end namespace m4d
