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
// -------------------------------------------------------------------------------- 
/*
    m4dGeodesicBS.h

  Copyright (c) 2010-2014  Thomas Mueller


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

/*!  \class  m4d::GeodesicBS
     \brief  Calculate geodesics with the Bulirsch-Stoer method.



 */
// -------------------------------------------------------------------------------- 

#ifndef M4D_GEODESIC_BS_H
#define M4D_GEODESIC_BS_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dGeodesic.h>

#include "m4dGeodesic.h"

#define  DEF_BS_TINY         1e-30
#define  DEF_BS_MAX_ROW_NUM  8
#define  DEF_BS_MAX_ROW_NUMS (DEF_BS_MAX_ROW_NUM+1)
#define  DEF_BS_SAFE1        0.25
#define  DEF_BS_SAFE2        0.7
#define  DEF_BS_MAX_RED      1e-5
#define  DEF_BS_MIN_RED      0.7
#define  DEF_BS_MAX_SCALE    0.1

namespace m4d
{

// ---------------------------------------------------
//    class definition:   GeodesicBS
// ---------------------------------------------------
class MOTION_API GeodesicBS : public Geodesic
{
public:
    GeodesicBS ( Metric* metric, enum_geodesic_type  type = enum_geodesic_lightlike );
    virtual ~GeodesicBS ();

    // --------- public methods -----------
public:
    virtual void setMaxAffineParamStep ( double hmax );
    void    setNumberOfSteps      ( int numSteps );

    virtual enum_break_condition  calculateGeodesic     ( const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &lambda );

    virtual enum_break_condition  calculateGeodesic     ( const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                                          vec4 *&points, vec4 *&dirs, int &numPoints );

    virtual enum_break_condition  calculateGeodesicData ( const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &epsilons, std::vector<double> &lambda );

    virtual enum_break_condition  calcParTransport      ( const vec4 initPos, const vec4 initDir,
                                                          const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                          const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs,
                                                          std::vector<double> &lambda,
                                                          std::vector<vec4> &base0, std::vector<vec4> &base1,
                                                          std::vector<vec4> &base2, std::vector<vec4> &base3 );

    virtual enum_break_condition  calcSachsJacobi       ( const vec4 initPos, const vec4 initCoordDir,
                                                          const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                                          const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
                                                          const enum_nat_tetrad_type  tetrad_type,
                                                          const int maxNumPoints,
                                                          std::vector<vec4> &points,  std::vector<vec4> &dirs,
                                                          std::vector<double> &lambda,
                                                          std::vector<vec4> &sachs0, std::vector<vec4> &sachs1,
                                                          std::vector<vec5> &jacobi, vec5 &maxJacobi );

    virtual enum_break_condition  calcSachsJacobi       ( const vec4 initPos, const vec4 initCoordDir,
                                                          const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                                          const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                          const enum_nat_tetrad_type  tetrad_type,
                                                          const int maxNumPoints,
                                                          vec4 *&points, vec4 *&dirs,
                                                          double *&lambda,
                                                          vec4 *&sachs0, vec4 *&sachs1,
                                                          vec5 *&jacobi, vec5 &maxJacobi, int &numPoints );


    virtual enum_break_condition  nextStep              ( double htry, double &hdid, double &hnext, double &constraint );
    virtual enum_break_condition  nextStepSachsJacobi   ( double htry, double &hdid, double &hnext, double &constraint );

    virtual bool    nextStep            ( int &status );
    virtual bool    nextStepPar         ( int &status );
    virtual bool    nextStepSachsJacobi ( int &status );


    virtual void                  print ( FILE* fptr = stderr );


    // --------- protected methods -----------
protected:
    void       modMidPoint    (  double* yy, double* dydx, double H, int numSteps, double* yout );
    void       modMidPointJS  (  double* yy, double* dydx, double H, int numSteps, double* yout );
    void       polyExtrpol    ( int iest, double xest, double* yest, double* yz, double* dy );

    // -------- protected attribute ---------
protected:
    double     ym[DEF_MAX_YS];
    double     yn[DEF_MAX_YS];
    double     dydx[DEF_MAX_YS];
    double     yerr[DEF_MAX_YS];
    double     ysav[DEF_MAX_YS];
    double     yseq[DEF_MAX_YS];
    double     yscal[DEF_MAX_YS];
    double     cc[DEF_MAX_YS];

    double**   dd;
    double     a[DEF_BS_MAX_ROW_NUMS];
    double     alf[DEF_BS_MAX_ROW_NUM][DEF_BS_MAX_ROW_NUM];
    double     err[DEF_BS_MAX_ROW_NUM];
    double     x[DEF_BS_MAX_ROW_NUM];
    double     eps,epsold,xnew;

    int        nSeq[DEF_BS_MAX_ROW_NUMS+1];

    int        mFirst, mKmax, mKopt;
    int        mReduct, mExitFlag;

    long double     mhmin;
    double     mhmax;
    int        mNumCoords;
    int        mNumSteps;
};

} // end namespace m4d

#endif
