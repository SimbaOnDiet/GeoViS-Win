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
    m4dFermiWalker.h  
 
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

/*!  \class  m4d::FermiWalker
     \brief  The local tetrad of an arbitrarily moving object has to be Fermi-Walker transported.
  
        The Fermi-Walker transport equation for the local tetrad \f$\mathbf{e}_{(i)}=e_{(i)}^{\mu}\partial_{\mu}\f$ is given by:
         \f[ \frac{de_{(i)}^{\mu}}{d\tau} = -\Gamma_{\alpha\beta}^{\mu}e_{(0)}^{\alpha}e_{(i)}^{\beta} + \frac{1}{c}\left(\eta_{(i)(k)}a^{(k)}e_{(0)}^{\mu}-\eta_{(0)(i)}a^{(k)}e_{(k)}^{\mu}\right) \f]


        There are two possibilities to calculate the Fermi-Walker transport of
        a local tetrad. On the one hand, the proper 3-acceleration with respect 
        to the local tetrad can be given for each time step. On the other hand,
        the worldline can be defined as a function \f$x=x(\tau)\f$. 

        Because the four-acceleration is always orthogonal to the four-velocity,
        the zero- or time- component of the proper acceleration with respect to
        the local tetrad where the e0- base vector is tangential to the four-velocity
        vanishes.

        The worldline function \f$x=x(\tau)\f$ has to take care of the 
        parameter pointer!


        Note that e0 is tangential to four-velocity : y[8]..y[11] = y[4]..y[7].
        
 */
// -------------------------------------------------------------------------------- 

#ifndef M4D_MOTION_FERMIWALKER_H
#define M4D_MOTION_FERMIWALKER_H

#include <iostream>

#include "m4dMotion.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   FermiWalker
// ---------------------------------------------------
class MOTION_API FermiWalker : public Motion
{
 public:
   FermiWalker ( Metric* metric );
   virtual ~FermiWalker ();
 
// --------- public methods -----------
 public:
   void  setCurrPropAccel   ( double a1, double a2, double a3 );
   void  getCurrPropAccel   ( double &a1, double &a2, double &a3 );
 
   bool  setInitialVelocity ( double fm, double v, double theta, double phi, 
                              enum_nat_tetrad_type  type = enum_nat_tetrad_default  );
   void  getInitialVelocity ( double &v, double &theta, double &phi );
   
   //! Reset proper time.
   void  resetProperTime    ( ) { mLambda = 0.0; }
 
   virtual enum_break_condition  calculateMotion( const vec4 initPos, double fm, double v, double theta_v, double phi_v,
                                                  double a, double theta_a, double phi_a, 
                                                  const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                  const int maxNumPoints,
                                                  std::vector<vec4> &points, 
                                                  std::vector<vec4> &base0, std::vector<vec4> &base1, std::vector<vec4> &base2, std::vector<vec4> &base3 );
 
 
   void  setCalcWithWorldline ( bool withworldline );

   void  set_x_tau       ( vec4 (*func)(double,void*) );
   void  set_u_tau       ( vec4 (*func)(double,void*) );
   void  set_a_tau       ( vec4 (*func)(double,void*) );
   void  set_params      ( void* params );

   bool  get_x_tau       ( double tau, vec4 &x );
   bool  get_u_tau       ( double tau, vec4 &u );
   bool  get_a_tau       ( double tau, vec4 &a );

   bool  initWorldline   ( double tauStart );
   bool  updateWorldline ( double tau );

   enum_break_condition  nextStep   ( );
   enum_break_condition  nextStepWL ( );


// --------- protected methods -----------
 protected:
   bool  calcDerivs   ( const double y[], double dydx[] );
   bool  calcDerivsWL ( const double y[], double dydx[] );
    

// -------- protected attribute ---------
 protected:
  //! Acceleration with respect to local tetrad.
  double mPropAcc[4];
  //! Current four-velocity in coordinates.
  double mInitVel[4];
  //! Current four-acceleration.
  vec4 mCurrAcc;
 
  double  mVel;
  double  mTheta;
  double  mPhi;

  bool     mCalcWithWorldline;
  vec4     (*x_tau)(double, void*);
  vec4     (*u_tau)(double, void*);
  vec4     (*a_tau)(double, void*);
  void*   mParams;
};

} // end namespace m4d

#endif

