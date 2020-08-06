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
    m4dObject.h

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

/*!  \class  m4d::Object
     \brief  Master object that stores all relevant data.


     Settings file consist of:
   \verbatim
        METRIC            <string>
        PARAM             <num>  <name>  <double>
        INIT_POS          <double>  <double>  <double>  <double>
        INIT_DIR          <double>  <double>  <double>
        INIT_ANGLE_VEL    <double>  <double>  <double>
        TIME_DIR          <int>
        GEOD_SOLVER_TYPE  <int>
        GEODESIC_TYPE     <int>
        STEPSIZE_CTRL     <int>
        STEPSIZE          <double>
        STEPSIZE_MAX_MIN  <double> <double>
        EPSILONS          <double> <double>
        CONSTR_EPSILON    <double>
        MAX_NUM_POINTS    <int>
        TETRAD_TYPE       <int>
        BASE_0            <double>  <double>  <double>  <double>
        BASE_1            <double>  <double>  <double>  <double>
        BASE_2            <double>  <double>  <double>  <double>
        BASE_3            <double>  <double>  <double>  <double>
        BOOST             <double>  <double>  <double>
        SPEED_OF_LIGHT    <double>
        GRAV_CONSTANT     <double>
        DIELECTRIC_PERM   <double>                                 \endverbatim

    Description:
    <ul>
      <li>METRIC   :  name of the metric defined in each class.
      <li>PARAM    :  continuous number,  name of the parameter,  value.
      <li>INIT_POS :  initial position in coordinates  (x0, x1, x2, x3).
      <li>INIT_DIR :  initial direction in coordinates (n1, n2, n3) with respect to local tetrad.
      <li>INIT_ANGLE_VEL : initial angles and velocity  (ksi, chi, vel) with respect to local tetrad.
      <li>TIME_DIR : time direction (>0: future, <0: past).
      <li>GEOD_SOLVER_TYPE : type of geodesic solver as numerical value.
      <li>GEODESIC_TYPE : type of geodesic (-1: timelike, 0: lightlike, 1: spacelike).
      <li>STEPSIZE_CTRL : stepsize control (1: yes, 0: no).
      <li>STEPSIZE : stepsize.
      <li>STEPSIZE_MAX_MIN : maximum and minimum stepsize.
      <li>EPSILONS : absolute and relative epsilons for geodesic integration.
      <li>CONSTR_EPSILON: epsilon for the constraint condition.
      <li>RESIZE_EPSILON: epsilon when resizing of y[4] shall be done.
      <li>MAX_NUM_POINTS : maximum number of points to be calculated.
      <li>TETRAD_TYPE : type of local tetrad.
      <li>BASE_0 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BASE_1 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BASE_2 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BASE_3 : base vector of local tetrad with respect to natural local tetrad at initial position.
      <li>BOOST  : boost parameters (ksi, chi, vel).
      <li>SPEED_OF_LIGHT : value of the speed of light (geometric units: 1, physical units: 299792458.0)
      <li>GRAV_CONSTANT  : value of the gravitational constant (geometric units: 1, physical units: 6.67...)
      <li>DIELECTRIC_PERM : value of the dielectric permittivity (geometric units: 1, physical units: ...)
    </ul>
    \sa testDatabase.cpp
    \sa enum_nat_tetrad_type
*/
// -------------------------------------------------------------------------------- 

#ifndef M4D_OBJECT_H
#define M4D_OBJECT_H


#include <iostream>
#include <fstream>
#include <cassert>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <motion/m4dMotionList.h>
#include <motion/m4dMotionDatabase.h>
#include <extra/m4dUtilities.h>

namespace m4d
{

// ---------------------------------------------------
//    class definition:   Object
// ---------------------------------------------------
class EXTRA_API Object
{
public:
    Object ( );
    ~Object ( );

    // --------- public methods -----------
public:
    void   clearAll           ( );
    void   resetAll           ( );

    bool   setLorentzTransf   ( const double chi, const double ksi, const double beta );
    void   resetLorentzTransf ( );

    bool   loadSettings       ( std::string filename, bool printset = false );
    bool   saveSettings       ( std::string filename, std::string dat = std::string());
    void   printSettings      ( FILE* fptr = stderr );

    bool   makeReport         ( std::string  &text );

    // --------- public attributes --------
public:
    MetricDatabase*       metricDB;
    Metric*               currMetric;
    IntegratorDatabase*   solverDB;
    Geodesic*             geodSolver;
    enum_integrator       geodSolverType;
    enum_geodesic_type    type;
    bool                  stepsizeControlled;
    double                stepsize;
    double                max_stepsize;
    double                min_stepsize;
    double                epsAbs;
    double                epsRel;
    double                epsConstr;
    double                epsResize;

    vec4                  startPos;
    vec3                  startDir;
    double                ksi;
    double                chi;
    double                vel;
    vec4                  coordDir;
    vec4                  base[4];
    bool                  isBaseInCoords;

    int                   axes_orient;
    double                boost_ksi;
    double                boost_chi;
    double                boost_beta;
    mat4                  lorentz;

    int                   timeDirection;
    enum_nat_tetrad_type  tetradType;
    unsigned int          maxNumPoints;
    std::vector<vec4>     points;
    std::vector<vec4>     dirs;
    std::vector<double>   lambda;
    std::vector<vec4>     sachs1;
    std::vector<vec4>     sachs2;
    std::vector<vec5>     jacobi;
    vec5                  maxJacobi;

    double                speed_of_light;
    double                grav_constant;
    double                dielectric_perm;
};

} // end namespace m4d

#endif
