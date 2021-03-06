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
   m4dMetricErnst.cpp

  Copyright (c) 2010-2014  Thomas Mueller


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

#include "m4dMetricErnst.h"

namespace m4d
{

#define eps 1.0e-6


/*! Standard constructor for the Ernst-Schwarzschild metric.
 *
 * \param  mass : mass of the black hole.
 * \param  B : magnetic field
 */
MetricErnst :: MetricErnst ( double mass, double B )
{
    mMetricName  = "Ernst";
    setCoordType(enum_coordinate_spherical);
    
    mPhysicalUnits = enum_physical_constants_geom;
    mSpeedOfLight   = 1.0;
    mGravConstant   = 1.0;
    mDielectricPerm = 1.0;
    
    addParam("mass",mass);
    mMass = mass;
    addParam("B",B);
    mB = B;

    mSign = 1.0;
    mLocTeds.push_back(enum_nat_tetrad_static);
    mDrawTypes.push_back(enum_draw_effpoti);

    setStandardValues();
}

/*!
 */
MetricErnst :: ~MetricErnst()
{
}


// *********************************** public methods ******************************

/*! Calculate the contravariant metric components at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricErnst :: calculateMetric  ( const double* pos )
{
    double r     = pos[1];
    double theta = pos[2];

    double c  = mSpeedOfLight;
    double rs = 2.0*mMass;
    double B  = mB;

    double t1 = c*c;
    double t3 = 1/r;
    double t5 = B*B;
    double t6 = r*r;
    double t7 = t5*t6;
    double t8 = sin(theta);
    double t9 = t8*t8;
    double t10 = t9*t1;
    double t17 = t5*t5;
    double t18 = t6*t6;
    double t19 = t17*t18;
    double t20 = t9*t9;
    double t21 = t20*t1;
    double t30 = 1/(1.0-rs*t3);
    double t47 = pow(1.0+t7*t9,2.0);

    g_compts[0][0] = -t1+t1*rs*t3-2.0*t7*t10+2.0*t5*r*t10*rs-t19*t21+t17*t6*r*t21*rs;
    g_compts[0][1] = 0.0;
    g_compts[0][2] = 0.0;
    g_compts[0][3] = 0.0;
    g_compts[1][0] = 0.0;
    g_compts[1][1] = t30+2.0*t7*t9*t30+t19*t20*t30;
    g_compts[1][2] = 0.0;
    g_compts[1][3] = 0.0;
    g_compts[2][0] = 0.0;
    g_compts[2][1] = 0.0;
    g_compts[2][2] = t6+2.0*t5*t18*t9+t17*t18*t6*t20;
    g_compts[2][3] = 0.0;
    g_compts[3][0] = 0.0;
    g_compts[3][1] = 0.0;
    g_compts[3][2] = 0.0;
    g_compts[3][3] = t6*t9/t47;

    return true;
}

/*! Calculate the Christoffel symbols of the second kind at position 'pos'.
 *
 *  \param pos : pointer to position.
 */
bool 
MetricErnst :: calculateChristoffels  ( const double* pos )
{
    double r     = pos[1];
    double theta = pos[2];

    double c  = mSpeedOfLight;
    double rs = 2.0*mMass;
    double B  = mB;

    double t1 = B*B;
    double t2 = r*r;
    double t3 = t2*r;
    double t5 = sin(theta);
    double t6 = t5*t5;
    double t8 = 4.0*t1*t3*t6;
    double t11 = t1*t6*rs*t2;
    double t13 = t8-3.0*t11+rs;
    double t14 = r-rs;
    double t16 = c*c;
    double t17 = t1*t2;
    double t18 = t17*t6;
    double t19 = 1.0+t18;
    double t20 = 1/t19;
    double t27 = 1/r;
    double t30 = cos(theta);
    double t35 = 1/t14;
    double t37 = t27*t20;
    double t39 = t13*t35*t37/2.0;
    double t40 = t5*t30;
    double t43 = 2.0*t17*t40*t20;
    double t55 = 3.0*t18+1.0;
    double t57 = t55*t27*t20;
    double t58 = -1.0+t18;
    double t60 = t58*t27*t20;
    double t63 = t58*t30;
    double t66 = t63*t20/t5;
    double t68 = t1*t1;
    double t69 = t2*t2;
    double t71 = t6*t6;
    double t74 = 1/(1.0+2.0*t18+t68*t69*t71);
    double t77 = t19*t19;
    double t79 = 1/t77/t19;

    christoffel[0][0][0] = 0.0;
    christoffel[0][0][1] = t13*t14*t16*t20/t3/2.0;
    christoffel[0][0][2] = 2.0*t14*t1*t27*t5*t16*t20*t30;
    christoffel[0][0][3] = 0.0;
    christoffel[0][1][0] = t39;
    christoffel[0][1][1] = 0.0;
    christoffel[0][1][2] = 0.0;
    christoffel[0][1][3] = 0.0;
    christoffel[0][2][0] = t43;
    christoffel[0][2][1] = 0.0;
    christoffel[0][2][2] = 0.0;
    christoffel[0][2][3] = 0.0;
    christoffel[0][3][0] = 0.0;
    christoffel[0][3][1] = 0.0;
    christoffel[0][3][2] = 0.0;
    christoffel[0][3][3] = 0.0;
    christoffel[1][0][0] = t39;
    christoffel[1][0][1] = 0.0;
    christoffel[1][0][2] = 0.0;
    christoffel[1][0][3] = 0.0;
    christoffel[1][1][0] = 0.0;
    christoffel[1][1][1] = (t8-5.0*t11-rs)*t35*t37/2.0;
    christoffel[1][1][2] = -2.0*t40*r*t1*t35*t20;
    christoffel[1][1][3] = 0.0;
    christoffel[1][2][0] = 0.0;
    christoffel[1][2][1] = t43;
    christoffel[1][2][2] = t57;
    christoffel[1][2][3] = 0.0;
    christoffel[1][3][0] = 0.0;
    christoffel[1][3][1] = 0.0;
    christoffel[1][3][2] = 0.0;
    christoffel[1][3][3] = -t60;
    christoffel[2][0][0] = t43;
    christoffel[2][0][1] = 0.0;
    christoffel[2][0][2] = 0.0;
    christoffel[2][0][3] = 0.0;
    christoffel[2][1][0] = 0.0;
    christoffel[2][1][1] = t43;
    christoffel[2][1][2] = t57;
    christoffel[2][1][3] = 0.0;
    christoffel[2][2][0] = 0.0;
    christoffel[2][2][1] = -t55*t14*t20;
    christoffel[2][2][2] = t43;
    christoffel[2][2][3] = 0.0;
    christoffel[2][3][0] = 0.0;
    christoffel[2][3][1] = 0.0;
    christoffel[2][3][2] = 0.0;
    christoffel[2][3][3] = -t66;
    christoffel[3][0][0] = 0.0;
    christoffel[3][0][1] = 0.0;
    christoffel[3][0][2] = 0.0;
    christoffel[3][0][3] = 0.0;
    christoffel[3][1][0] = 0.0;
    christoffel[3][1][1] = 0.0;
    christoffel[3][1][2] = 0.0;
    christoffel[3][1][3] = -t60;
    christoffel[3][2][0] = 0.0;
    christoffel[3][2][1] = 0.0;
    christoffel[3][2][2] = 0.0;
    christoffel[3][2][3] = -t66;
    christoffel[3][3][0] = 0.0;
    christoffel[3][3][1] = t14*t74*t6*t58*t79;
    christoffel[3][3][2] = t74*t5*t63*t79;
    christoffel[3][3][3] = 0.0;

    return true;
}

/*! Calculate Jacobi matrix. 
 *
 *  \param pos : pointer to position.
 */
bool 
MetricErnst :: calculateChrisD ( const double* pos )
{
    double r     = pos[1];
    double theta = pos[2];

    double c  = mSpeedOfLight;
    double rs = 2.0*mMass;
    double B  = mB;

    double t1 = c*c;
    double t2 = B*B;
    double t3 = r*r;
    double t4 = t3*t3;
    double t6 = sin(theta);
    double t7 = t6*t6;
    double t9 = 4.0*t2*t4*t7;
    double t10 = t2*t2;
    double t13 = t7*t7;
    double t15 = 4.0*t10*t4*t3*t13;
    double t16 = t3*r;
    double t17 = t2*t16;
    double t18 = t7*rs;
    double t19 = t17*t18;
    double t20 = 4.0*t19;
    double t23 = t13*rs;
    double t24 = t10*t4*r*t23;
    double t26 = t2*t7;
    double t27 = rs*rs;
    double t29 = t26*t27*t3;
    double t30 = 2.0*t29;
    double t31 = t10*t4;
    double t33 = t31*t13*t27;
    double t36 = 2.0*r*rs;
    double t40 = t2*t3;
    double t41 = t40*t7;
    double t43 = pow(1.0+t41,2.0);
    double t44 = 1/t43;
    double t49 = 1/r;
    double t52 = cos(theta);
    double t53 = r-rs;
    double t54 = t53*t53;
    double t63 = 2.0*t17*t7;
    double t65 = t26*rs*t3;
    double t69 = 1/t3;
    double t70 = t69*t44;
    double t76 = t52*t52;
    double t78 = t40*t7*t76;
    double t80 = -t76+t78+t7+t40*t13;
    double t81 = t80*t44;
    double t90 = 1/t54;
    double t93 = (-t9+t15+12.0*t19-6.0*t24-6.0*t29+3.0*t33+t36-t27)*t90*t70/2.0;
    double t98 = 4.0*t6*r*t2*t52*t44;
    double t100 = 2.0*t40*t81;
    double t114 = t2*r;
    double t120 = t31*t13;
    double t121 = 3.0*t120;
    double t124 = (t121+1.0)*t69*t44;
    double t128 = (-4.0*t41+t120-1.0)*t69*t44;
    double t130 = t114*t18;
    double t141 = t31*t13*t76;
    double t144 = t13*t7*t10*t4;
    double t148 = (-4.0*t78+t141-t7+t144-t76)*t44/t7;
    double t157 = 2.0*t41;
    double t159 = 1/(1.0+t157+t120);
    double t160 = t43*t43;
    double t161 = 1/t160;
    double t162 = t159*t161;

    chrisD[0][0][0][0] = 0.0;
    chrisD[0][0][0][1] = 0.0;
    chrisD[0][0][0][2] = 0.0;
    chrisD[0][0][0][3] = 0.0;
    chrisD[0][0][1][0] = 0.0;
    chrisD[0][0][1][1] = -t1*(-t9+t15+t20-14.0*t24-t30+9.0*t33+t36-3.0*t27)*t44/t4/2.0;
    chrisD[0][0][1][2] = 4.0*t2*t49*t6*t52*t54*t1*t44;
    chrisD[0][0][1][3] = 0.0;
    chrisD[0][0][2][0] = 0.0;
    chrisD[0][0][2][1] = -2.0*t2*t6*t1*t52*(t63-rs-3.0*t65)*t70;
    chrisD[0][0][2][2] = -2.0*t53*t2*t1*t81*t49;
    chrisD[0][0][2][3] = 0.0;
    chrisD[0][0][3][0] = 0.0;
    chrisD[0][0][3][1] = 0.0;
    chrisD[0][0][3][2] = 0.0;
    chrisD[0][0][3][3] = 0.0;
    chrisD[0][1][0][0] = 0.0;
    chrisD[0][1][0][1] = -t93;
    chrisD[0][1][0][2] = t98;
    chrisD[0][1][0][3] = 0.0;
    chrisD[0][1][1][0] = 0.0;
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
    chrisD[0][2][0][1] = t98;
    chrisD[0][2][0][2] = -t100;
    chrisD[0][2][0][3] = 0.0;
    chrisD[0][2][1][0] = 0.0;
    chrisD[0][2][1][1] = 0.0;
    chrisD[0][2][1][2] = 0.0;
    chrisD[0][2][1][3] = 0.0;
    chrisD[0][2][2][0] = 0.0;
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
    chrisD[0][3][3][0] = 0.0;
    chrisD[0][3][3][1] = 0.0;
    chrisD[0][3][3][2] = 0.0;
    chrisD[0][3][3][3] = 0.0;
    chrisD[1][0][0][0] = 0.0;
    chrisD[1][0][0][1] = -t93;
    chrisD[1][0][0][2] = t98;
    chrisD[1][0][0][3] = 0.0;
    chrisD[1][0][1][0] = 0.0;
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
    chrisD[1][1][0][0] = 0.0;
    chrisD[1][1][0][1] = 0.0;
    chrisD[1][1][0][2] = 0.0;
    chrisD[1][1][0][3] = 0.0;
    chrisD[1][1][1][0] = 0.0;
    chrisD[1][1][1][1] = -(-t9+t15+t20-10.0*t24-t30+5.0*t33-t36+t27)*t90*t70/2.0;
    chrisD[1][1][1][2] = t98;
    chrisD[1][1][1][3] = 0.0;
    chrisD[1][1][2][0] = 0.0;
    chrisD[1][1][2][1] = 2.0*t52*t6*t2*(t63+rs-t65)*t90*t44;
    chrisD[1][1][2][2] = 2.0*t114*t80/t53*t44;
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
    chrisD[1][2][1][1] = t98;
    chrisD[1][2][1][2] = -t100;
    chrisD[1][2][1][3] = 0.0;
    chrisD[1][2][2][0] = 0.0;
    chrisD[1][2][2][1] = -t124;
    chrisD[1][2][2][2] = t98;
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
    chrisD[1][3][3][1] = t128;
    chrisD[1][3][3][2] = -t98;
    chrisD[1][3][3][3] = 0.0;
    chrisD[2][0][0][0] = 0.0;
    chrisD[2][0][0][1] = t98;
    chrisD[2][0][0][2] = -t100;
    chrisD[2][0][0][3] = 0.0;
    chrisD[2][0][1][0] = 0.0;
    chrisD[2][0][1][1] = 0.0;
    chrisD[2][0][1][2] = 0.0;
    chrisD[2][0][1][3] = 0.0;
    chrisD[2][0][2][0] = 0.0;
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
    chrisD[2][1][1][1] = t98;
    chrisD[2][1][1][2] = -t100;
    chrisD[2][1][1][3] = 0.0;
    chrisD[2][1][2][0] = 0.0;
    chrisD[2][1][2][1] = -t124;
    chrisD[2][1][2][2] = t98;
    chrisD[2][1][2][3] = 0.0;
    chrisD[2][1][3][0] = 0.0;
    chrisD[2][1][3][1] = 0.0;
    chrisD[2][1][3][2] = 0.0;
    chrisD[2][1][3][3] = 0.0;
    chrisD[2][2][0][0] = 0.0;
    chrisD[2][2][0][1] = 0.0;
    chrisD[2][2][0][2] = 0.0;
    chrisD[2][2][0][3] = 0.0;
    chrisD[2][2][1][0] = 0.0;
    chrisD[2][2][1][1] = -(8.0*t41+t121-4.0*t130+1.0)*t44;
    chrisD[2][2][1][2] = -4.0*t40*t6*t52*t53*t44;
    chrisD[2][2][1][3] = 0.0;
    chrisD[2][2][2][0] = 0.0;
    chrisD[2][2][2][1] = t98;
    chrisD[2][2][2][2] = -t100;
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
    chrisD[2][3][3][1] = -t98;
    chrisD[2][3][3][2] = t148;
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
    chrisD[3][0][3][0] = 0.0;
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
    chrisD[3][1][3][1] = t128;
    chrisD[3][1][3][2] = -t98;
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
    chrisD[3][2][3][1] = -t98;
    chrisD[3][2][3][2] = t148;
    chrisD[3][2][3][3] = 0.0;
    chrisD[3][3][0][0] = 0.0;
    chrisD[3][3][0][1] = 0.0;
    chrisD[3][3][0][2] = 0.0;
    chrisD[3][3][0][3] = 0.0;
    chrisD[3][3][1][0] = 0.0;
    chrisD[3][3][1][1] = -(7.0*t120-8.0*t10*t16*t23-12.0*t41+12.0*t130+1.0)*t7*t162;
    chrisD[3][3][1][2] = -2.0*(t121-6.0*t41+1.0)*t53*t6*t52*t159*t161;
    chrisD[3][3][1][3] = 0.0;
    chrisD[3][3][2][0] = 0.0;
    chrisD[3][3][2][1] = -4.0*(t157-3.0)*t7*t6*t52*t114*t162;
    chrisD[3][3][2][2] = -(t144+7.0*t141-12.0*t78+t76-t7)*t159*t161;
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
MetricErnst :: localToCoord ( const double* pos, const double* ldir, double* dir,
                              enum_nat_tetrad_type   )
{ 
    double r     = pos[1];
    double theta = pos[2];

    double w = sqrt(1.0 - 2.0*mMass/r);
    double L = 1.0 + mB*mB*r*r*sin(theta)*sin(theta);

    dir[0] = ldir[0]/L/w;
    dir[1] = ldir[1]*w/L;
    dir[2] = ldir[2]/(L*r);
    dir[3] = ldir[3]*L/(r*sin(theta));
}

/*! Transform coordinate 4-direction to local 4-direction.
 *
 *  \param  pos  :  pointer to position array.
 *  \param  cdir :  pointer to coordinate direction.
 *  \param  ldir :  pointer to calculated local direction array.
 *  \param  type :  type of tetrad.
 */
void 
MetricErnst :: coordToLocal ( const double* pos, const double* cdir, double* ldir,
                              enum_nat_tetrad_type   )
{
    double r     = pos[1];
    double theta = pos[2];

    double w = sqrt(1.0 - 2.0*mMass/r);
    double L = 1.0 + mB*mB*r*r*sin(theta)*sin(theta);

    ldir[0] = cdir[0]*L*w;
    ldir[1] = cdir[1]*L/w;
    ldir[2] = cdir[2]*L*r;
    ldir[3] = cdir[3]*r*sin(theta)/L;
}


/*! Test break condition.
 *
 *  \param pos    : pointer to position array.
 *  \return true  : radial position r < 0.0 or  r^2<=(1.0+eps)*rs^2.
 *  \return false : position is valid.
 */
bool 
MetricErnst :: breakCondition ( const double* pos)
{
    bool br = false;
    if ((pos[1]<0.0) || (pos[1]*pos[1]<=(1.0+eps)*4.0*mMass*mMass)) { br=true; }
    
    return br;
}


/*! Tests whether the constraint equation is fulfilled.
 *
 *  double the constraint equation for lightlike and timelike geodesics reads:
 \verbatim
     sum = g_{\mu\nu} dot(x)^{\mu} dot(x)^{\nu} - kappa c^2 = 0.
 \endverbatim
 *  \param  y[]   : pointer to position and direction coordinates.
 *  \param  kappa : timelike (-1.0), lightlike (0.0).
 *  \return double : sum.
 */
double 
MetricErnst :: testConstraint ( const double y[], const double kappa )
{
    double r     = y[1];
    double theta = y[2];

    double w2 = 1.0 - 2.0*mMass/r;
    double L  = 1.0 + mB*mB*r*r*sin(theta)*sin(theta);

    // Scale the directions with the speed of light before doubling them !!
    double dt     = y[4];
    double dr     = y[5];
    double dtheta = y[6];
    double dphi   = y[7];

    double sum = -kappa*mSign;
    sum += L*L*(-w2*dt*dt + dr*dr/w2 + r*r*dtheta*dtheta) + r*r*sin(theta)*sin(theta)/(L*L)*dphi*dphi;
    return sum;
}

/*! Set parameter 'pName' to 'val'.
 *
 *  Set 'mass' parameter and adjust Schwarzschild radius  rs=2GM/c^2.
 *  'charge' represents the charge of the black hole.
 */
bool 
MetricErnst :: setParam ( std::string pName, double val )
{
    Metric::setParam(pName,val);

    if (pName=="mass")
    {
        mMass = val;
    }
    else if (pName=="b")
        mB = val;
    
    return true;
}

/*! Effective potential.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param type : geodesic type.
 *  \param x : abscissa value.
 *  \param val : reference to effective potential value.
 *  \return true : effective potential exists at x.
 */
bool
MetricErnst :: effPotentialValue  ( const vec4 pos, const vec4 cdir ,enum_geodesic_type type, const double x, double &val )
{
    double kappa = 0.0;
    if (type==enum_geodesic_timelike)
        kappa = mSign;

    if (pos[1]<2.0*mMass+1e-2 || x<2.0*mMass+1e-2)
        return false;

    double Lc = 1.0+mB*mB*pos[1]*pos[1]*sin(pos[2])*sin(pos[2]);
    double h  = pow(pos[1]*sin(pos[2])/Lc,2.0)*cdir[3];
    double k  = (1.0-2.0*mMass/pos[1])*Lc*Lc*cdir[0];

    double Lambda = 1.0 + mB*mB*x*x;
    val = h*h*(1.0-2.0*mMass/x)/(x*x) - k*k*pow(Lambda,-4.0) + kappa*(1.0-2.0*mMass/x)*pow(Lambda,-2.0);
    return true;
}


/*! Total energy.
 *  \param pos : initial position.
 *  \param cdir : initial four-direction.
 *  \param x : abscissa value.
 *  \param val : reference to total energy value.
 *  \return true : effective potential exists at x.
 */
bool
MetricErnst :: totEnergy ( const vec4 , const vec4 , const double , double &val )
{
    val = 0.0;
    return true;
}

/*! Generate report.
 */
bool    
MetricErnst :: report ( const vec4 pos, const vec4 cdir, std::string &text )
{
    std::stringstream ss;
    ss << "Report for the  Ernst metric\n\tcoordinate : (t,r,theta,phi)\n";
    ss << "---------------------------------------------------------------\n";
    ss.precision(DEF_FIXED_REPORT_PRECISION);
    ss.setf(std::ios::fixed);

    double Lc = 1.0+mB*mB*pos[1]*pos[1]*sin(pos[2])*sin(pos[2]);
    double h  = pow(pos[1]*sin(pos[2])/Lc,2.0)*cdir[3];
    double k  = (1.0-2.0*mMass/pos[1])*Lc*Lc*cdir[0];
    ss << "  constant of motion ........... h = " << h << std::endl;
    ss << "  constant of motion ........... k = " << k << std::endl;

    text = ss.str();
    return true;
}

// ********************************* protected methods *****************************
/*!
 */
void 
MetricErnst :: setStandardValues( )
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

} // end namespace m4d
