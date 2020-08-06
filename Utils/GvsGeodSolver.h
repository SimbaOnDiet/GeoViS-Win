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

#ifndef GVS_GEOD_SOLVER_H
#define GVS_GEOD_SOLVER_H

#include <iostream>
#include <cstdio>

#include "Obj/GvsBase.h"
#include "Obj/STMotion/GvsLocalTetrad.h"
#include <motion/m4dGeodesic.h>
#include <motion/m4dMotionList.h>
#include <metric/m4dMetric.h>

class GvsGeodSolver : public GvsBase
{
public:
    GvsGeodSolver( m4d::Metric* metric,
                   m4d::enum_integrator m4dGeodSolver );

    GvsGeodSolver( m4d::Metric* metric, m4d::enum_geodesic_type type, m4d::enum_time_direction tDir,
                   m4d::enum_integrator m4dGeodSolver );

    ~GvsGeodSolver();

    void setMetric( m4d::Metric* metric );
    m4d::Metric* getMetric();

    bool setSolver( m4d::enum_integrator m4dGeodSolver );

    void                     setGeodType( m4d::enum_geodesic_type gType );
    m4d::enum_geodesic_type  getGeodType() const;

    void                      setTimeDir ( m4d::enum_time_direction tDir );
    m4d::enum_time_direction  getTimeDir ( void  ) const;

    void   setEpsilons( double eps_a, double eps_r );
    void   getEpsilons( double &eps_a, double &eps_r );

    void   setStepSizeControl( const bool cn );
    bool   getStepSizeControl();

    void   setStepsize( double step );
    double getStepsize() const;

    void   setMaxStepsize( double step );
    double getMaxStepsize() const;


    int startConditionLocal ( const m4d::vec4* pos, m4d::vec4 &dir );

    int startCondition ( const m4d::vec4* pos, m4d::vec4 &dir );
    int startCondition ( const double pos[], double dir[] );

    double calcNormCondition ( const double pos[], const double dir[] );
    double calcNormCondition ( const m4d::vec4 pos, const m4d::vec4 dir );


    m4d::enum_break_condition calculateGeodesic ( const m4d::vec4& yStart, const m4d::vec4& Dir,
                                                  const double maxNumPoints,
                                                  m4d::vec4 *&points, m4d::vec4 *&dirs, int &numPoints );

    m4d::enum_break_condition calcParTransport ( const m4d::vec4& yStart, const m4d::vec4& yDir, const m4d::vec4 base[4],
                                                 const double maxNumPoints,
                                                 GvsLocalTetrad *&lt, int &numPoints );

    m4d::enum_break_condition calcSachsJacobi ( const m4d::vec4 &startOrig, const m4d::vec4 &startDir, const m4d::vec3 &localDir,
                                                const int maxNumPoints,
                                                const GvsLocalTetrad *lt,
                                                m4d::vec4 *&points, m4d::vec4 *&dirs, double *&lambda,
                                                m4d::vec4 *&sachs1, m4d::vec4 *&sachs2,
                                                m4d::vec5 *&rayJacobi, m4d::vec5 &rayMaxJacobi, int &numPoints );


    void  setBoundingBox  ( const double p1[4], const double p2[4] );
    void  setBoundingBox  ( const m4d::vec4 p1, const m4d::vec4 p2 );
    void  getBoundingBox  ( double p1[4], double p2[4] );
    void  setBoundingTime ( double minTime,  double maxTime );
    void  getBoundingTime ( double &minTime, double &maxTime );

    void errorMessage ( m4d::enum_break_condition brCond ) const;

    void Print ( FILE* fptr = stderr );

protected:
    bool   outsideBoundingBox ( const double* pos );

private:
    m4d::Metric*     mMetric;
    m4d::enum_geodesic_type  mGeodType;
    m4d::enum_time_direction mTimeDir;

    m4d::Geodesic* m4dSolver;
    m4d::enum_integrator m4dGeodSolverType;
    std::string solverName;

    double epsilon_abs;
    double epsilon_rel;

    bool   stepSizeControlled;
    double stepSize;       //!< Initial stepsize for calculate Geodesic
    double maxStepsize;    //!< Maximum stepsize for calculation


    // The bounding box parameters for the four coordinates
    double  boundBoxMin[4];
    double  boundBoxMax[4];

};

#endif
