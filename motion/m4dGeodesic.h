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
    m4dGeodesic.h

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

/*!  \class  m4d::Geodesic
     \brief  Base class for all geodesics.

       The Geodesic class cannot be used itself due to the pure virtual
       methods calculateGeodesic() and calculateGeodesicData().

       The geodesic equation reads
       \f[ \frac{d^2x^{\mu}}{d\lambda^2} + \Gamma_{\alpha\beta}^{\mu}\frac{dx^{\alpha}}{d\lambda}\frac{dx^{\beta}}{d\lambda} = 0, \f]
       where \f$\lambda\f$ is an affine parameter and \f$\Gamma_{\alpha\beta}^{\mu}\f$ are the Christoffell symbols of the second
       kind defined by
       \f[ \Gamma_{\alpha\beta}^{\mu} = \frac{1}{2}g^{\mu\rho}\left(g_{\rho\alpha,\beta}+g_{\rho\beta,\alpha}-g_{\alpha\beta,\rho}\right). \f]
       In the case of a timelike geodesic, this affine parameter equals the proper time.
       A timelike \f$(\kappa=-1)\f$ or a lightlike \f$(\kappa=0)\f$ geodesic has to fulfill the constraint equation
       \f[ g_{\mu\nu}\frac{dx^{\mu}}{d\lambda}\frac{dx^{\nu}}{d\lambda}=\kappa c^2.\f]

       The parallel transport equation for a vector X along the geodesic with tangent u reads
       \f[ \frac{dX^{\mu}}{d\lambda} + \Gamma_{\alpha\beta}^{\mu}u^{\alpha}X^{\beta} = 0. \f]
*/
// -------------------------------------------------------------------------------- 

#ifndef M4D_GEODESIC_H
#define M4D_GEODESIC_H


#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include <m4dGlobalDefs.h>
#include <metric/m4dMetric.h>
#include <motion/m4dMotion.h>

namespace m4d
{

// ---------------------------------------------------
//    class definition:   Geodesic
// ---------------------------------------------------
class MOTION_API Geodesic : public Motion
{
public:
    Geodesic ( Metric* metric, enum_geodesic_type  type = enum_geodesic_lightlike );
    virtual ~Geodesic ();

    // --------- public methods -----------
public:
    void    setGeodesicType           ( enum_geodesic_type  type );
    enum_geodesic_type   type         ( );

    void    setEpsilons               ( double eps_a, double eps_r);
    void    getEpsilons               ( double &eps_a, double &eps_r );

    void    setStepSizeControlled     ( bool control = true );

    void    setCalcWithParTransport   ( bool calcwith = false );
    bool    calcWithParTransport      ( );

    void    setResize                 ( double eps, double factor );
    void    getResize                 ( double &eps, double &factor );

    virtual enum_break_condition  initializeGeodesic    ( const vec4 initPos, const vec4 initDir, double &cstr );

    virtual enum_break_condition  calculateGeodesic     ( const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &lambda ) = 0;

    virtual enum_break_condition  calculateGeodesic     ( const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                                          vec4 *&points, vec4 *&dirs, int &numPoints ) = 0;

    virtual enum_break_condition  calculateGeodesicData ( const vec4 initPos, const vec4 initDir, const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs, std::vector<double> &epsilons, std::vector<double> &lambda ) = 0;

    virtual enum_break_condition  calcParTransport      ( const vec4 initPos, const vec4 initDir,
                                                          const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                          const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs,
                                                          std::vector<double> &lambda,
                                                          std::vector<vec4> &base0, std::vector<vec4> &base1,
                                                          std::vector<vec4> &base2, std::vector<vec4> &base3 ) = 0;

    virtual enum_break_condition  calcSachsJacobi       ( const vec4 initPos, const vec4 initCoordDir,
                                                          const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                                          const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                          const enum_nat_tetrad_type  tetrad_type,
                                                          const int maxNumPoints,
                                                          std::vector<vec4> &points, std::vector<vec4> &dirs,
                                                          std::vector<double> &lambda,
                                                          std::vector<vec4> &sachs0, std::vector<vec4> &sachs1,
                                                          std::vector<vec5> &jacobi, vec5 &maxJacobi ) = 0;

    virtual enum_break_condition  calcSachsJacobi       ( const vec4 initPos, const vec4 initCoordDir,
                                                          const vec3 localNullDir, const vec3 locX, const vec3 locY, const vec3 locZ,
                                                          const vec4 e0, const vec4 e1, const vec4 e2, const vec4 e3,
                                                          const enum_nat_tetrad_type  tetrad_type,
                                                          const int maxNumPoints,
                                                          vec4 *&points, vec4 *&dirs,
                                                          double *&lambda,
                                                          vec4 *&sachs0, vec4 *&sachs1,
                                                          vec5 *&jacobi, vec5 &maxJacobi, int &numPoints ) = 0;

    virtual bool    nextStep            ( int &status ) = 0;
    virtual bool    nextStepPar         ( int &status ) = 0;
    virtual bool    nextStepSachsJacobi ( int &status ) = 0;

    virtual double  testConstraint      ( );

    virtual void    print ( FILE* fptr = stderr );

    // --------- protected methods -----------
protected:
    void  setKappa              ( );
    virtual bool  calcDerivs            ( const double y[], double dydx[] );
    virtual bool  calcDerivsPar         ( const double y[], double dydx[] );
    virtual bool  calcDerivsSachsJacobi ( const double y[], double dydx[] );

    void  calcSachsBasis        ( const vec3 localNullDir, const vec3 locX, const vec3 loxY, const vec3 locZ );
    void  setSachsBasis         ( const vec4 s1, const vec4 s2 );
    void  calcJacobiParams      ( const double lambda, const double y[], vec5 &currJacobi );

    void  findMaxJacobi         ( vec5 &currJacobi, vec5 &maxJacobi );

    // -------- protected attribute ---------
protected:
    //! Type of geodesic.
    enum_geodesic_type  mType;
    double              mKappa;

    bool       mCalcWithParTransport;
    int        mNumCoords;

    bool       mStepsizeControlled;
    //! Absolute epsilon.
    double     epsilon_abs;
    //! Relative epsilon.
    double     epsilon_rel;
    //! Epsilon for resize of y[4].
    double     resizeEps;
    double     resizeFac;

    vec3       mSachsBasisB1;
    vec3       mSachsBasisB2;
};

} // end namespace m4d

#endif
