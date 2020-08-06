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
    m4dMetricPravda_C_Can.h

  Copyright (c) 2010-2014  Thomas Mueller, Frank Grave, Felix Beslmeisl


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

/*!  \class  m4d::MetricPravda_C_Can
     \brief  This is the C-Metric as given in class MetricPTD_C in "canonical coordinates adapted to the boost-rotation symmetry".

             The line element is given by

             \f[ds^2 = \frac{1}{z^2-t^2}\left(e^{\rho}r^2(z\,dt-t\,dz)^2 -e^{\lambda}(z\,dz-t\,dt)^2 \right) - e^{\lambda}\,dr^2-r^2e^{-\rho}\,d\varphi^2\f]
             with \f[e^{\rho}=\frac{R_3+R+Z_3-r^2}{4\alpha^2\left(R_1+R+Z_1-r^2\right)},\f]
                  \f[e^{\lambda}=\frac{2\alpha^2\left(R(R+R_1+Z_1)-Z_1r^2\right) \left(R_1R_3+(R+Z_1)(R+Z_3)-(Z_1+Z_3)r^2) \right)}{R_iR_3\left( R(R +R_3+Z_3) - Z_3r^2\right)},\f]
                  \f[R=\frac12\left(z^2-t^2+r^2\right),\f]
                  \f[R_i=\sqrt{(R + Z_i)^2 - 2Z_i r^2},\f]
                  \f[Z_i=z_i-z_2,\f]
                  \f[\alpha^2=\frac14 \frac{m^2}{A^6(z_2-z_1)^2(z3-z1)^2},\f]
                  \f[q=\frac1{4\alpha^2},\f]
             and \f$z_3<z_1<z2\f$ the roots of \f$2A^4z^3-A^2z^2+m^2.\f$
             A and m are real parameters.

             The natural local tetrad for \f$z^2-t^2>0\f$ is given by
             \f[  \mathbf{e}_{(t)} = \frac{1}{\sqrt{z^2-t^2}}\left( qze^{-\rho/2}\,\partial_t, + te^{-\lambda/2}\,\partial_z,\right)\quad
                  \mathbf{e}_{(r)} = e^{-\lambda/2}\,\partial_r, \quad
                  \mathbf{e}_{(\varphi)} = re^{\rho/2}\,\partial_{\varphi}, \quad
                  \mathbf{e}_{(z)} = \frac{1}{\sqrt{z^2-t^2}}\left( qte^{-\rho/2}\,\partial_t, + ze^{-\lambda/2}\,\partial_z,\right)\f]

             and for \f$z^2-t^2<0\f$ by
             \f[  \mathbf{e}_{(t)} = \frac{1}{\sqrt{z^2-t^2}}\left( qte^{-\rho/2}\,\partial_t, + ze^{-\lambda/2}\,\partial_z,\right) \quad
                  \mathbf{e}_{(r)} = e^{-\lambda/2}\,\partial_r, \quad
                  \mathbf{e}_{(\varphi)} = re^{\rho/2}\,\partial_{\varphi}, \quad
                  \mathbf{e}_{(z)} = \frac{1}{\sqrt{z^2-t^2}}\left( qze^{-\rho/2}\,\partial_t, + te^{-\lambda/2}\,\partial_z,\right).\f]

             This metric ist dicussed in<br>
             Pravda, V. and Pravdov&aacute; ,A., <b>Co-accelerated particles in the C-metric</b>  Class. Quantum Gravit. <b>18</b>, 1205 (2001).
*/

// --------------------------------------------------------------------------------
#ifndef M4DMETRICPRAVDA_C_CAN_H
#define M4DMETRICPRAVDA_C_CAN_H

#include "m4dMetric.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   MetricPravda_C_Can
// ---------------------------------------------------
class MetricPravda_C_Can : public Metric
{
 public:
  MetricPravda_C_Can ( double A = 0.01,  double m = 1.0);
  virtual ~MetricPravda_C_Can ();

// --------- public methods -----------
 public:
  virtual bool   calculateMetric        ( const double* pos );
  virtual bool   calculateChristoffels  ( const double* pos );
  virtual bool   calculateChrisD        ( const double* pos );

  virtual void   localToCoord   ( const double* pos, const double* ldir, double* dir,
                                  enum_nat_tetrad_type  type = enum_nat_tetrad_default );
  virtual void   coordToLocal   ( const double* pos, const double* cdir, double* ldir,
                                  enum_nat_tetrad_type  type = enum_nat_tetrad_default );

  virtual bool   breakCondition ( const double* pos);

  virtual double testConstraint ( const double y[], const double kappa );

  virtual bool   setParam ( std::string pName, double val );

  virtual void   calculateRoots(vec3 & roots, double p, double q);

  virtual bool   report ( const vec4 pos, const vec4 cdir, std::string &text );


// --------- specific public methods ----------
 public:

  double rho(double tau, double zeta, double eta);
  double rho_tau(double tau, double zeta, double eta);
  double rho_zeta(double tau, double zeta, double eta);
  double rho_eta(double tau, double zeta, double eta);

  double lambda(double tau, double zeta, double eta);
  double lambda_tau(double tau, double zeta, double eta);  //diff(lambda(tau,zeta,eta),tau)
  double lambda_zeta(double tau, double zeta, double eta); //diff(lambda(tau,zeta,eta),zeta)
  double lambda_eta(double tau, double zeta, double eta);  //diff(lambda(tau,zeta,eta),eta)

// --------- protected methods -----------
 protected:
  virtual void setStandardValues ( );


// -------- protected attribute ---------
 protected:
  double Par_A;
  double Par_m;
  vec3 z_i;   //z_i Roots of Polynom
  double Z1; //z_1-z_2
  double Z3; //z_3-2_2
  double alpha2;
  double q;

};

} // end namespace m4d

#endif //
