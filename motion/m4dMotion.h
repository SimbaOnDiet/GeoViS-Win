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
    m4dMotion.h  
 
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

/*!  \class  m4d::Motion
     \brief  Base class for motion/geodesics in 4d spacetimes.
  
             The Motion class cannot be used itself due to the pure virtual
             method calcDerivs().

             The initial position as well as the initial direction have to
             be given in coordinates of the corresponding metric. For the 
             initial direction, the localToCoord() method of the corresponding
             Metric class can be used. The initial position is checked whether it
             violates the break condition of the metric.

             For internal calculations, the current position, direction, and
             base vectors are mapped to the double array y:

            \verbatim
               y[ 0...3]  :  current position,
               y[ 4...7]  :  current direction,
               y[ 8..11]  :  e0,
               y[12..15]  :  e1,
               y[16..19]  :  e2.
               y[20..23]  :  e3. \endverbatim
             
        
*/
// -------------------------------------------------------------------------------- 

#ifndef M4D_MOTION_H
#define M4D_MOTION_H

#include <iostream>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <extra/m4dUtilities.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

namespace m4d
{
  
// ---------------------------------------------------
//    class definition:   Motion
// ---------------------------------------------------
class METRIC_API Motion
{
 public:
   Motion ( Metric* metric );
   virtual ~Motion ();
 
// --------- public methods -----------
 public:
   virtual bool  setInitialPosition  ( double ipos[4] );
   virtual bool  setInitialPosition  ( vec4 ipos );
 
   virtual bool  setInitialDirection ( double d0, double d1, double d2, double d3 );
   virtual bool  setInitialDirection ( vec4 dir );

   virtual void  setInitialTetrad    ( double e0[4], double e1[4], double e2[4], double e3[4] );
   virtual void  setInitialTetrad    ( vec4 e0, vec4 e1, vec4 e2, vec4 e3 );

   virtual bool  setLocalInitialDir  ( double d1, double d2, double d3, double beta = 1.0 );
   virtual bool  setLocalInitialDir  ( vec3 dir, double beta = 1.0 );

   virtual  void        setMetric    ( Metric* metric );
   virtual  Metric*     getMetric    ( );

   void    getPosition           ( double* p );
   vec4    getPosition           ( );

   void    getDirection          ( double* p );
   vec4    getDirection          ( );
 
   void    resetAffineParam      ( double pt = 0.0 );
   double  getAffineParam        ( );
   void    setAffineParamStep    ( double step );
   double  getAffineParamStep    ( );
   void    resetAffineParamStep  ( );

   virtual void    setMaxAffineParamStep ( double step );
           double  getMaxAffineParamStep ( );
   virtual void    setMinAffineParamStep ( double step );
           double  getMinAffineParamStep ( );

   void    setE                  ( unsigned int i, vec4 ee );
   vec4    getE                  ( unsigned int i );
   void    getE0                 ( double* p );
   void    getE1                 ( double* p );
   void    getE2                 ( double* p );
   void    getE3                 ( double* p );
   void    getTetrad             ( double* e0, double* e1, double* e2, double* e3 );
   void    getTetrad             ( vec4 &e0, vec4 &e1, vec4 &e2, vec4 &e3 );
   void    getTetrad             ( mat4 &m );
   void    getTetradInv          ( mat4 &m );
   void    getTetrad             ( float* m );
   void    getTetradInv          ( float* m );
        
   vec4    coordToLocal          ( vec4 cv );

   bool    isOrthonormal         ( );

   void    setBoundingBox        ( double p1[4], double p2[4] );
   void    setBoundingBox        ( vec4 p1, vec4 p2 );
   void    getBoundingBox        ( vec4 &p1, vec4 &p2 );
   bool    outsideBoundBox       ( );
   
   void    setConstrEps          ( double eps );
   double  getConstrEps          ( );
  
   bool    isRightHanded         ( );
   bool    gramSchmidtOrth       ( );

   virtual double  testConstraint ( );

   void    getCurrentArray       ( double* cy );

   //! Get duration of calculation.
   double  getCalcTime ( ) { return mCalcTime; }
    
   void    printTetrad ( FILE* fptr = stderr );

// --------- protected methods -----------
 protected:
   //! Calculate the right side of the parallel transport.
   virtual bool calcDerivs       ( const double yn[], double dydx[] ) = 0;
 
// -------- protected attribute ---------
 protected:
   //! Pointer to the actual metric.
   Metric*    mMetric;

   //! Affine parameter.
   double     mLambda;
   //! Affine parameter stepsize.
   double     mLambdaStep;
   //! Initial affine parameter stepsize.
   double     mLambdaStepInit;
   //! Maximum affine parameter stepsize.
   double     mMaxLambdaStep;
   //! Minimum affine parameter stepsize.
   double     mMinLambdaStep;

   //! This array holds the current position, direction, and all tetrad vectors.
   double     y[DEF_MAX_YS];  
   //! Epsilon for the constraint equation.
   double     mConstraintEpsilon;
   
   //! Bounding box minimum.
   double     mBoundBoxMin[4];
   //! Bounding box maximum.
   double     mBoundBoxMax[4];

   //! Time in seconds for calculation of geodesic.
   double     mCalcTime;
};

} // end namespace m4d

#endif

