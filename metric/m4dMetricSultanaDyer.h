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
    m4dMetricSultanaDyer.h
 
  Copyright (c) 2010-2014  Thomas Mueller
 
 
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

/*!  \class  m4d::MetricSultanaDyer
     \brief  SultanaDyer cosmological black hole in spherical coordinates (t,r,theta,phi).

      The line element is given by
      \f[ds^2 = t^4\left[\left(1-\frac{2m}{r}\right)dt^2-\frac{4m}{r}dt\,dr-\left(1+\frac{2m}{r}\right)dr^2-r^2\left(d\vartheta^2+\sin^2\vartheta\,d\varphi^2\right)\right].\f]
      
      The comoving natural local tetrad reads
      \f[ \mathbf{e}_{(0)} = \frac{1+2m/r}{t^2\sqrt{1+2m/r}}\partial_t - \frac{2m/r}{t^2\sqrt{1+2m/r}}\partial_r,\quad \mathbf{e}_{(1)} = \frac{1}{t^2\sqrt{1+2m/r}}\partial_r,\quad \mathbf{e}_{(2)} = \frac{1}{t^2r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{t^2r\sin\vartheta}\partial_{\varphi}.\f]

      The static natural local tetrad reads
      \f[ \mathbf{e}_{(0)} = \frac{1}{t^2\sqrt{1-2m/r}},\quad \mathbf{e}_{(1)} = \frac{2m/r}{t^2\sqrt{1-2m/r}}\partial_t + \frac{\sqrt{1-2m/r}}{t^2}\partial_r, \quad \mathbf{e}_{(2)} = \frac{1}{t^2r}\partial_{\vartheta},\quad \mathbf{e}_{(3)} = \frac{1}{t^2r\sin\vartheta}\partial_{\varphi}.\f]

      For further information see...<br>
      J. Sultana and C.C. Dyer, "Cosmological black holes: A black hole in the Einstein-de Sitter universe", Gen. Relativ. Gravit. <b>37</b>, 1349--1370 (2005).
 */
// -------------------------------------------------------------------------------

#ifndef M4D_METRIC_SULTANA_DYER_COS_BH_H
#define M4D_METRIC_SULTANA_DYER_COS_BH_H

#include "m4dMetric.h"

namespace m4d
{

// ---------------------------------------------------
//    class definition:   MetricSultanaDyer
// ---------------------------------------------------
class MetricSultanaDyer: public Metric
{
  public:
    MetricSultanaDyer ( double mass = 1.0 );
    ~MetricSultanaDyer ();


// --------- public methods --------------
  public:
    virtual bool   calculateMetric       ( const double* pos );
    virtual bool   calculateChristoffels ( const double* pos );
    virtual bool   calculateChrisD       ( const double* pos );

    virtual void   localToCoord          ( const double* pos, const double* ldir, double* dir,
                                           enum_nat_tetrad_type type =enum_nat_tetrad_default );
    virtual void   coordToLocal          ( const double* pos, const double* cdir, double* ldir,
                                           enum_nat_tetrad_type type = enum_nat_tetrad_default );

    virtual bool   breakCondition        ( const double* pos);

    virtual double testConstraint        ( const double y[], const double kappa );

    virtual bool   setParam  ( std::string pName, double val );

    virtual bool   report    ( const vec4 pos, const vec4 cdir, std::string &text );

// --------- protected methods -----------
  protected:
    virtual void setStandardValues ();


// -------- protected attribute ----------
  protected:
    double mMass;
    
};

} // end namespace m4d

#endif
