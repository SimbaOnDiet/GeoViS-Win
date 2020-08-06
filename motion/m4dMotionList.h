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
    m4dMotionList.h

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
// --------------------------------------------------------------------------------
#ifndef  M4D_MOTION_LIST_H
#define  M4D_MOTION_LIST_H

#include <motion/m4dMotion.h>
#include <motion/m4dGeodesic.h>
#include <motion/m4dGeodesicBS.h>
#include <motion/m4dGeodesicGSL.h>
#include <motion/m4dGeodesicRK4.h>
#include <motion/m4dGeodesicDP54.h>
#include <motion/m4dGeodesicDP65.h>

namespace m4d
{

#ifdef USE_DP_INT
const int NUM_GEOD_SOLVERS = 12;
#else
const int NUM_GEOD_SOLVERS = 10;
#endif
const int GSL_INT_OFFSET   = 4;


/* ----------------------------------------------
 *   List of all solvers currently implemented
 *
 *   When editing this list please take care 
 *   of the ordering.
 * ---------------------------------------------- */
static const char  stl_solver_names[NUM_GEOD_SOLVERS][60] = {
  "Runge-Kutta (fourth order)",
  "Bulirsch-Stoer",
  "GSL: Embedded 2nd order Runge-Kutta",
  "GSL: 4th order (classical) Runge-Kutta",
  "GSL: Embedded 4th order Runge-Kutta-Fehlberg",
  "GSL: Embedded 4th order Runge-Kutta Cash-Karp",
  "GSL: Embedded 8th order Runge-Kutta Prince-Dormand",
  "GSL: Implicit 2nd order Runge-Kutta at Gaussian points",
  //"GSL: Implicit Burlisch-Stoer of Bader and Deuflhard",
  "GSL: M=1 implicit Gear",
  "GSL: M=2 implicit Gear"
#ifdef USE_DP_INT
  ,"Dormand-Prince: RK5(4)",
  "Dormand-Prince: RK6(5)"
#endif
};

static const char  stl_solver_nicknames[NUM_GEOD_SOLVERS][30] = {
  "RK4",
  "BS",
  "GSL_RK2",
  "GSL_RK4",
  "GSL_RK_Fehlberg",
  "GSL_RK_Cash-Karp",
  "GSL_RK_Prince-Dormand",
  "GSL_RK_Gaussian-points",
 // "GSL_BS_Bader-Deuflhard",
  "GSL_M1_IG",
  "GSL_M2_IG"
#ifdef USE_DP_INT
  ,"DP54",
  "DP65"
#endif
};

enum enum_integrator
{
  gsUnknown = -1,
  gsIrk4 = 0,       // standard fourth-order Runge-Kutta
  gsInrbs,          // Bulirsch-Stoer
  gsIgslrk2,        // gsl_odeiv_step_rk2
  gsIgslrk4,        // gsl_odeiv_step_rk4
  gsIgslfehlberg,   // gsl_odeiv_step_rkf45
  gsIgslcash,       // gsl_odeiv_step_rkck
  gsIgslprinc,      // gsl_odeiv_step_rk8pd
  gsIi2,            // gsl_odeiv_step_rk2imp
//  gsIbs,            // gsl_odeiv_step_bsimp
  gsIm1,            // gsl_odeiv_step_gear1
  gsIm2             // gsl_odeiv_step_gear2
#ifdef USE_DP_INT
 ,gsIdp54,          // Dormand-Prince: RK5(4)
  gsIdp65           // Dormand-Prince: RK6(5)
#endif
};

} // end namespace m4d

#endif

