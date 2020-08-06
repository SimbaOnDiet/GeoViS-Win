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
    m4dGeodesicGSL.h

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

/*!  \class  m4d::GeodesicGSL
     \brief  Calculate geodesics with the GNU Scientific Library

             GSL:

 */
// -------------------------------------------------------------------------------- 

#ifndef M4D_GEODESIC_GSL_H
#define M4D_GEODESIC_GSL_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "m4dGeodesic.h"

#define  DEF_GSL_LAMBDA_MAX  1.0e38

//extern const gsl_odeiv_step_type*  bugfix_step_type;

namespace m4d
{

enum enum_gslint_type	
{
    enum_gslint_geodesic = 0,
    enum_gslint_geodesic_data,
    enum_gslint_partrans,
    enum_gslint_partrans_jacobi
};

// prototpyes
int func_adaptor_geod   ( double x, const double y[], double f[], void *params );
int jac_adaptor_geod    ( double x, const double y[], double *dfdy, double dfdt[], void *params );

int func_adaptor_par    ( double x, const double y[], double f[], void *params );

int func_adaptor_jacobi ( double x, const double y[], double f[], void *params );


// ---------------------------------------------------
//    class definition:   GeodesicGSL
// ---------------------------------------------------
class MOTION_API GeodesicGSL : public Geodesic
{
public:
    GeodesicGSL ( Metric* metric, const gsl_odeiv_step_type*  step_type, int solver_type,
                  enum_geodesic_type  type = enum_geodesic_lightlike );
    virtual ~GeodesicGSL ();

    // --------- public methods -----------
public:
    virtual enum_break_condition  initializeGeodesic    ( const vec4 initPos, const vec4 initDir, double &cstr );

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

    virtual enum_break_condition  calcSachsJacobi       (const vec4 initPos, const vec4 initCoordDir,
                                                          const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                                          const vec4 b0, const vec4 b1, const vec4 b2, const vec4 b3,
                                                          const enum_nat_tetrad_type  tetrad_type,
                                                          const int maxNumPoints,
                                                          vec4 *&points, vec4 *&dirs,
                                                          double *&lambda,
                                                          vec4 *&sachs0, vec4 *&sachs1,
                                                          vec5 *&jacobi, vec5 &maxJacobi, int &numPoints );

    int   func_geod   ( double x, const double y[], double f[], void *params );
    int   jac_geod    ( double x, const double y[], double *dfdy, double dfdt[], void *params );

    int   func_par    ( double x, const double y[], double f[], void *params );

    int   func_jacobi ( double x, const double y[], double f[], void *params );

    virtual bool  nextStep            ( int &status );
    virtual bool  nextStepPar         ( int &status );
    virtual bool  nextStepSachsJacobi ( int &status );

    virtual void  print ( FILE* fptr = stderr );


    // --------- protected methods -----------
protected:
    void initialize  ( const gsl_odeiv_step_type*  step_type,  enum_gslint_type type = enum_gslint_geodesic );
    void allocMemory ( );
    void freeMemory  ( );

    // -------- protected attribute ---------
protected:
    const gsl_odeiv_step_type*  mStepType;
    int                         mSolverType;
    gsl_odeiv_step*      mStep;
    gsl_odeiv_control*   mControl;
    gsl_odeiv_evolve*    mEvolve;

    gsl_odeiv_system     mSys;

    double    dydt_in[DEF_MAX_YS], dydt_out[DEF_MAX_YS];
};

} // end namespace m4d

#endif

