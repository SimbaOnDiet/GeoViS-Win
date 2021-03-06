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
   m4dMetricErezRosenVar.cpp

  Copyright (c) 2010-2014   Thomas Mueller


   double this file is part of the m4d-library.

   double the m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   double the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   double the m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dMetricErezRosenVar.h"

namespace m4d
{

MetricErezRosenVar :: MetricErezRosenVar ( double mass, double q )
{
    mMetricName  = "ErezRosenVar";
    setCoordType(enum_coordinate_spherical);
    
    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight   = 1.0;
    mGravConstant   = 1.0;
    mDielectricPerm = 1.0;
    
    addParam("mass",mass);
    mMass = mass;
    addParam("q",q);
    mQ = q;

    mSign = 1.0;
    mLocTeds.push_back(enum_nat_tetrad_static);
    setStandardValues();
}


MetricErezRosenVar :: ~MetricErezRosenVar() {
}


// *********************************** public methods ******************************


bool 
MetricErezRosenVar :: calculateMetric  ( const double* pos )
{
    double r     = pos[1];
    double theta = pos[2];

    double m = mMass;
    double psi,g,Delta;
    if (!calcPotentials(pos,psi,g,Delta)) {
        return false;
    }

    double t1 = psi;         // psi(r,theta);
    double t2 = exp(t1);
    double t3 = t2*t2;
    double t4 = g;           // g(r,theta);
    double t5 = exp(t4);
    double t6 = t5*t5;
    double t7 = 1/t3;
    double t8 = t6*t7;
    double t9 = Delta;       // Delta(r,theta);
    double t10 = r*r;
    double t11 = m*r;
    double t18 = sin(theta);
    double t19 = t18*t18;
    double t20 = t7*t19;

    g_compts[0][0] = -t3;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t8*t9/(t10-2.0*t11);
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t8*t9;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t20*t10-2.0*t20*t11;

    return true;
}


bool 
MetricErezRosenVar :: calculateChristoffels  ( const double* pos )
{
    double r     = pos[1];
    double theta = pos[2];

    double m = mMass;
    double psi,g,Delta;
    double dpsidr, dpsidtheta, dgdr,dgdtheta, dDdr, dDdtheta;
    calcDiffPots(pos,psi,g,Delta,dpsidr, dpsidtheta, dgdr,dgdtheta, dDdr, dDdtheta);

    double t2 = -r+2.0*m;
    double t3 = r*t2;
    double t4 = psi;              // psi(r,theta);
    double t5 = exp(t4);
    double t6 = t5*t5;
    double t7 = t6*t6;
    double t9 = g;                // g(r,theta);
    double t10 = exp(t9);
    double t11 = t10*t10;
    double t12 = 1/t11;
    double t13 = Delta;           // Delta(r,theta);
    double t14 = 1/t13;
    double t15 = t12*t14;
    double t16 = dpsidr;          // diff(psi(r,theta),r);
    double t20 = dpsidtheta;      // diff(psi(r,theta),theta);
    double t25 = 1/r/t2;
    double t26 = dgdr;            // diff(g(r,theta),r);
    double t27 = t13*t26;
    double t28 = r*r;
    double t31 = m*r;
    double t34 = t13*t16;
    double t39 = dDdr;            // diff(Delta(r,theta),r);
    double t52 = dgdtheta;        // diff(g(r,theta),theta);
    double t57 = dDdtheta;        // diff(Delta(r,theta),theta);
    double t59 = t14*(-2.0*t13*t52+2.0*t13*t20-t57);
    double t62 = t59/2.0;
    double t66 = t14*(-2.0*t27+2.0*t34-t39);
    double t67 = t66/2.0;
    double t72 = -t16*t28+r+2.0*t16*r*m-m;
    double t73 = t25*t72;
    double t76 = sin(theta);
    double t78 = cos(theta);
    double t81 = (t76*t20-t78)/t76;
    double t83 = t76*t76;
    double t87 = t76*r;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = -t3*t7*t15*t16;
    christoffel[0][0][2] = t7*t12*t14*t20;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t16;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t20;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t16;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = -t25*t14*(2.0*t27*t28-4.0*t27*t31-2.0*t34*t28+4.0*t34*t31+t39*t28-2.0*t39*r*m-2.0*t13*r+2.0*t13*m)/2.0;
    christoffel[1][1][2] = -t59*t25/2.0;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = -t62;
    christoffel[1][2][2] = -t67;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t73;
    christoffel[2][0][0] = t20;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = -t62;
    christoffel[2][1][2] = -t67;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t3*t66/2.0;
    christoffel[2][2][2] = -t62;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = -t81;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t73;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = -t81;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t3*t12*t14*t83*t72;
    christoffel[3][3][2] = -t15*t87*(-t87*t20+r*t78+2.0*t76*m*t20-2.0*m*t78);
    christoffel[3][3][3] = 0.0;
    return true;
}


bool 
MetricErezRosenVar :: calculateChrisD ( const double* )
{
    /*
    double r     = pos[1];
    double theta = pos[2];

    // TODO   
    */
    return true;
}


void
MetricErezRosenVar :: localToCoord ( const double* pos, const double* ldir, double* dir,
                                     enum_nat_tetrad_type   )
{ 
    double r     = pos[1];
    double theta = pos[2];

    double m = mMass;
    double psi,g,Delta;
    calcPotentials(pos,psi,g,Delta);

    dir[0] = ldir[0]*exp(-psi);
    dir[1] = ldir[1]*exp(psi-g)*sqrt(r*r-2*m*r)/sqrt(Delta);
    dir[2] = ldir[2]*exp(psi-g)/sqrt(Delta);
    dir[3] = ldir[3]*exp(psi)/(sqrt(r*r-2*m*r)*sin(theta));
}


void 
MetricErezRosenVar :: coordToLocal ( const double* , const double* , double* ,
                                     enum_nat_tetrad_type   )
{
    //double r     = pos[1];
    //double theta = pos[2];

    // TODO
    //ldir[0] = cdir[0]*L*w;
    //ldir[1] = cdir[1]*L/w;
    //ldir[2] = cdir[2]*L*r;
    //ldir[3] = cdir[3]*r*sin(theta)/L;
}


bool 
MetricErezRosenVar :: breakCondition ( const double* )
{
    bool br = false;
    // if ((pos[1]<0.0) || (pos[1]*pos[1]<=(1.0+eps)*4.0*mMass*mMass)) { br=true; }
    // TOOO
    return br;
}


double 
MetricErezRosenVar :: testConstraint ( const double y[], const double kappa )
{
    double r     = y[1];
    double theta = y[2];

    double m = mMass;
    double psi,g,Delta;
    calcPotentials(y,psi,g,Delta);

    double dt     = y[4];
    double dr     = y[5];
    double dtheta = y[6];
    double dphi   = y[7];

    double st = sin(theta);

    double sum = -kappa*mSign;
    sum += -exp(2*psi)*dt*dt + exp(2*g-2*psi)*Delta*(dr*dr/(r*r-2*m*r)+dtheta*dtheta) + exp(-2*psi)*(r*r-2*m*r)*st*st*dphi*dphi;
    return sum;
}


bool 
MetricErezRosenVar :: setParam ( std::string pName, double val )
{
    Metric::setParam(pName,val);
    if (pName=="mass") {
        mMass = val;
    }
    else if (pName=="q") {
        mQ = val;
    }
    return true;
}


bool    
MetricErezRosenVar :: report ( const vec4 pos, const vec4 , std::string &text )
{
    std::stringstream ss;
    ss << "Report for the  Erez-Rosen metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);


    ss <<"Velocity for timelike circular orbits: beta = sqrt[f/(1+f)]" << std::endl;
    ss <<"with f=(m(-2m^3q+m^2qr-5r^3))/(-8m^4q+6m^3q(2m-r)+4m^3qr+15mr^3-5r^4)"<< std::endl;
    ss <<std::endl;
    ss <<"Mass:                                     m="<<mMass<<std::endl;
    ss <<"Position:                                 r="<<pos[1]<<std::endl;
    ss <<"Charge:                                   q="<<mQ<<std::endl;
    ss <<"Velocity for timelike circular orbits: beta="<<getCircularVelocity(pos[1])<<std::endl;

    // TODO
    text = ss.str();
    return true;
}


/*! Determine the velocity for a closed circular orbit if it exists.
 *   A circular timelike geodesic with respect to r-coordinate does exist
 *   only for r>=3rs (last timelike circular orbit).
 * \param r  Radial coordinate.
 * \param tedType type of tetrad.
 */
double
MetricErezRosenVar :: getCircularVelocity( const double r, const enum_nat_tetrad_type  )
{

    double m = mMass;
    double q = mQ;

    double m2=m*m;
    double m3=m*m*m;
    double m4=m*m*m*m;

    double r3=r*r*r;
    double r4=r*r*r*r;

    double f=(m*(-2.0*m3*q+m2*q*r-5.0*r3))/(-8.0*m4*q+6.0*m3*q*(2.0*m-r)+4.0*m3*q*r+15.0*m*r3-5.0*r4);

    return sqrt(f/(1.0+f));
}

vec4
MetricErezRosenVar :: getCircularFourVel( const vec4 pos, const enum_nat_tetrad_type )
{
    double beta = getCircularVelocity(pos[1]);
    if (beta>0.0 && beta<1.0) {
        double gamma = 1.0/sqrt(1.0-beta*beta);
        vec4 e0,e1,e2,e3;
        getNatTetrad(pos,e0,e1,e2,e3);
        return mSpeedOfLight * gamma*(e0 + beta*e3);
    }
    return vec4();
}

// ********************************* protected methods *****************************

void 
MetricErezRosenVar :: setStandardValues( )
{
    mInitPos[0] = 0.0;
    mInitPos[1] = 6.0;
    mInitPos[2] = M_PI_2;
    mInitPos[3] = 0.0;
    mInitDir[0] = 1.0;
    mInitDir[1] = 0.0;
    mInitDir[2] = 0.0;

    mCoordNames[0] = std::string("t");
    mCoordNames[1] = std::string("r");
    mCoordNames[2] = std::string("theta");
    mCoordNames[3] = std::string("phi");
}

bool
MetricErezRosenVar::calcPotentials( const double* pos, double &psi, double &g, double &delta ) {
    double r = pos[1];
    double theta = pos[2];
    double m = mMass;
    //double q = mQ;

    double st = sin(theta);
    delta = r*r - 2.0*m*r + m*m*st*st;
    g     = 0.5*log((r*r-2.0*m*r)/delta);
#if 1
    psi = 0.5*log(1.0-2.0*m/r);
#else
    double P2 = 0.5*(3*cos(theta)*cos(theta)-1.0);
    psi = 0.5*log(1.0-2.0*m/r) - 2.0/15.0*q*pow(m/r,3.0)*P2;
#endif
    return true;
}

bool
MetricErezRosenVar::calcDiffPots( const double* pos, double &psi, double &g, double &delta,
                                  double &dpsidr, double &dpsidtheta,
                                  double &dgdr, double &dgdtheta,
                                  double &dDdr, double &dDdtheta) {
    double r = pos[1];
    double theta = pos[2];
    double m = mMass;
    //double q = mQ;

    calcPotentials(pos,psi,g,delta);

    double st = sin(theta);
    double ct = cos(theta);

    dDdr = 2.0*r - 2.0*m;
    dDdtheta = 2.0*m*m*st*ct;

    dgdr = ((m*m*r-m*m*m)*st*st)/((m*m*r*r-2*m*m*m*r)*st*st+r*r*r*r-4*m*r*r*r+4*m*m*r*r);
    dgdtheta = -(m*m*ct*st)/(m*m*st*st+r*r-2*m*r);

#if 1
    dpsidr = m/(r*r-2.0*m*r);
    dpsidtheta = 0.0;
#else
    dpsidr = (m*m*m*q*(3*ct*ct-1))/(5*r*r*r*r)+m/((1-(2*m)/r)*r*r);
    dpsidtheta = (2*m*m*m*q*ct*st)/(5*r*r*r);
#endif
    return true;
}

} // end namespace m4d
