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
    m4dGeodesicDP65.h
 
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

/*!  \class  m4d::GeodesicDP65
     \brief  Calculate geodesics with a Dormand-Prince [RK6(5)] method.
  
         The algorithm is taken from
         Andreas Guthmann,
         "Einfuehrung in die Himmelsmechanik und Ephemeridenrechnung",
         Spektrum Verlag (2000), 2. Auflage, Seite 277

         Order of the method is p=5;
        
 */
// -------------------------------------------------------------------------------- 

#ifndef M4D_GEODESIC_DP65_H
#define M4D_GEODESIC_DP65_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dGeodesic.h>

#ifdef USE_DP_INT

namespace m4d
{

// ---------------------------------------------------
//    class definition:   GeodesicDP65
// ---------------------------------------------------
class MOTION_API GeodesicDP65 : public Geodesic
{
 public:
   GeodesicDP65 ( Metric* metric, enum_geodesic_type  type = enum_geodesic_lightlike );
   virtual ~GeodesicDP65 ();
 
// --------- public methods -----------
 public:
   void    setStepSizeControlled ( bool control = true );

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
                                                         const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                         const enum_nat_tetrad_type  tetrad_type,
                                                         const int maxNumPoints,
                                                         std::vector<vec4> &points, std::vector<vec4> &dirs, 
                                                         std::vector<double> &lambda, 
                                                         std::vector<vec4> &sachs0, std::vector<vec4> &sachs1, 
                                                         std::vector<vec5> &jacobi, vec5 &maxJacobi );

    virtual enum_break_condition  calcSachsJacobi       ( const vec4 initPos, const vec4 initCoordDir,
                                                          const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                                          const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
                                                          const enum_nat_tetrad_type  tetrad_type,
                                                          const int maxNumPoints,
                                                          vec4 *&points, vec4 *&dirs,
                                                          double *&lambda,
                                                          vec4 *&sachs0, vec4 *&sachs1,
                                                          vec5 *&jacobi, vec5 &maxJacobi, int &numPoints );

   //void  nextStep              ( double* yo, double* yn, double* yerr, double h );
   //void  nextStepPar           ( double* yo, double* yn, double* yerr, double h );
   //void  nextStepSachsJacobi   ( double* yo, double* yn, double* yerr, double h );
   
   virtual bool    nextStep            ( int &status );
   virtual bool    nextStepPar         ( int &status );
   virtual bool    nextStepSachsJacobi ( int &status );

   virtual void                  print ( FILE* fptr = stderr );   
   
// -------- protected attribute ---------
 protected: 
   bool       mCalcWithParTransport;
   int        mNumCoords;

   double a21,a31,a32,a41,a42,a43,a51,a52,a53,a54;
   double a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a82,a83,a84,a85,a86,a87;
   double b1,b2,b3,b4,b5,b6,b7,b8,db1,db2,db3,db4,db5,db6,db7,db8;

   int    mOrder;
   double stepSigma;
   double stepFac;

   double yn[DEF_MAX_YS],yerr[DEF_MAX_YS],h;
};

} // end namespace m4d

#endif // USE_DP_INT
#endif

