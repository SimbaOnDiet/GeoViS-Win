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
// -------------------------------------------------------------------------------
/*
    m4dMetricSchwarzschild.h

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

/*!  \class  m4d::MetricSchwarzschild
     \brief  Schwarzschild metric in spherical coordinates (t,r,theta,phi).

             The line element is given by

             \f[ds^2 = -\left(1-\frac{r_s}{r}\right) c^2dt^2 + \frac{dr^2}{1-r_s/r} + r^2\left(d\vartheta^2 + \sin^2\vartheta\,d\varphi^2\right),\f]
             where \f$r_s=2GM/c^2\f$ is the Schwarzschild radius. G is Newton's
             constant, M is the mass of the black hole, and c is the speed of light.

             The natural local tetrad is given by
             \f[ \mathbf{e}_{(0)} = \frac{1}{c\,\sqrt{1-r_s/r}}\partial_t,\quad \mathbf{e}_{(1)}=\sqrt{1-\frac{r_s}{r}}\partial_r,\quad \mathbf{e}_{(2)}=\frac{1}{r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{r\sin\vartheta}\partial_{\varphi}.\f]

             The embedding diagram (Flamm's paraboloid) is defined by the function
             \f[ z = 2\sqrt{r_s}\sqrt{r-r_s}. \f]

             From the Euler-Lagrange formalism, we have
             \f[ \frac{1}{2}\dot{r}^2+V_{\mbox{eff}} = \frac{1}{2}\frac{k^2}{c^2},\qquad\mbox{with}\qquad V_{\mbox{eff}}=\frac{1}{2}\left(\frac{h^2}{r^2}-\kappa c^2\right)\left(1-\frac{r_s}{r}\right)\f]
             and the constants of motion \f$k=(1-r_s/r)c^2\dot{t}\f$ and \f$h=r^2\dot{\varphi}\f$.

             Detailed discussions about the Schwarzschild metric can be found
             in the standard literature, e.g. \ref lit_wald "Wald", \ref lit_rindler "Rindler", \ref lit_mtw "MTW".

*/
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SCHWARZSCHILD_H
#define M4D_METRIC_SCHWARZSCHILD_H

#include "m4dMetric.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   MetricSchwarzschild
// ---------------------------------------------------
class MetricSchwarzschild : public Metric
{
public:
    MetricSchwarzschild ( double mass = 1.0 );
    virtual ~MetricSchwarzschild ();

    // --------- public methods -----------
public:
    virtual bool   calculateMetric        ( const double* pos );
    virtual bool   calculateChristoffels  ( const double* pos );
    virtual bool   calculateChrisD        ( const double* pos );

    virtual bool   calculateRiemann       ( const double* pos );
    virtual bool   calculateWeyl          ( const double* pos );
    virtual bool   calculateRicRotCoeffs  ( const double* pos );
    virtual bool   calculateContrRRC      ( const double* pos );

    virtual void   localToCoord           ( const double* pos, const double* ldir, double* dir,
                                            enum_nat_tetrad_type  type = enum_nat_tetrad_default );
    virtual void   coordToLocal           ( const double* pos, const double* cdir, double* ldir,
                                            enum_nat_tetrad_type  type = enum_nat_tetrad_default );

    virtual bool   breakCondition         ( const double* pos);

    virtual bool   calcDerivs             ( const double y[], double dydx[] );
    virtual bool   calcDerivsPar          ( const double y[], double dydx[] );
    virtual bool   calcDerivsSachsJacobi  ( const double y[], double dydx[] );
    virtual bool   calcDerivsFW           ( const double a[], const double y[], double dydx[] );

    virtual double testConstraint         ( const double y[], const double kappa );

    virtual bool   calcProduct            ( const double* pos, const double* u, const double* v, double &prod, bool preCalcMetric = true );

    virtual bool   setParam               ( std::string pName, double val );

    virtual bool   transToEmbedding       ( vec4 p, vec4 &ep );
    virtual bool   transToTwoPlusOne      ( vec4 p, vec4 &cp );
    virtual bool   transToCustom          ( vec4 p, vec4 &cp );

    virtual bool   setEmbeddingParam      ( std::string name, double val );
    virtual bool   testEmbeddingParams    ( );
    virtual int    getEmbeddingVertices   ( std::vector<vec3> &verts,
                                            std::vector<int> &indices, unsigned int &numElems, unsigned int &counter );

    virtual bool   effPotentialValue      ( const vec4 pos, const vec4 cdir, enum_geodesic_type type, const double x, double &val );
    virtual bool   totEnergy              ( const vec4 pos, const vec4 cdir, const double x, double &val );

    virtual double getCircularVelocity    ( const double r, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default );
    virtual vec4   getCircularFourVel     ( const vec4 pos, const enum_nat_tetrad_type  tedType = enum_nat_tetrad_default );

    virtual void   usePhysicalUnits       ( const enum_physical_constants  units );
    virtual void   setUnits               ( const double speed_of_light, const double grav_const, const double diel_perm );

    virtual bool   report ( const vec4 pos, const vec4 cdir, std::string &text );


    // --------- specific public methods ----------
public:
    bool  calcKsiCrit ( const vec4 pos, double &ksicrit );


    // --------- protected methods -----------
protected:
    virtual void  setStandardValues  ( );
    virtual void  contrChrisVecVec   ( const double y[], const double v[], const double w[], double* z, bool calc = true );
    virtual void  contrChrDVecVecVec ( const double y[], const double u[], const double v[], const double w[], double* z, bool calc = true );

    // -------- protected attribute ---------
protected:
    //! Schwarzschild radius rs=2GM/c^2.
    double rs;
    double mMass;

    double mEmb_rmin;
    double mEmb_rmax;
    double mEmb_rstep;
    double mEmb_phistep;
    double mEmb_r_num;
    double mEmb_phi_num;
};

} // end namespace m4d

#endif
