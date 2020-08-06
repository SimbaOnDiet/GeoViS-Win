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
//READING DONE
// ---------------------------------------------------------------------
//  Copyright (c) 2013-2014, Universitaet Stuttgart, VISUS, Thomas Mueller
//
//  This file is part of GeoViS.
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
// ---------------------------------------------------------------------
#ifndef GVS_BOUND_BOX_H
#define GVS_BOUND_BOX_H

#include <iostream>
#include <float.h>

#include "GvsGlobalDefs.h"

class GvsBoundBox
{
  public:
    GvsBoundBox ( void );
    GvsBoundBox ( const m4d::vec3& lowBounds, const m4d::vec3& uppBounds, bool testBounds = true );
    GvsBoundBox ( double x,  double y,  double z,
                  double dx, double dy, double dz, bool testBounds = true);
    GvsBoundBox ( const GvsBoundBox& otherBox );

    virtual ~GvsBoundBox(){}

    bool      isEmpty   ( ) const;
    m4d::vec3       lowBounds ( ) const;
    m4d::vec3       uppBounds ( ) const;
    m4d::vec3       size      ( ) const;
    double    size      ( int coord ) const;   // Tiefe, Hoehe
    double    volume    ( ) const;
    double    surface   ( ) const;

    bool      contains  ( const m4d::vec3 &point ) const;
    void      extendBoxToContain ( const m4d::vec3& pt );

    // -------- operators --------
    GvsBoundBox&  operator=  ( const GvsBoundBox &otherBox );
    GvsBoundBox   operator+  ( const GvsBoundBox &otherBox ) const;   // union
    GvsBoundBox&  operator+= ( const GvsBoundBox &otherBox );
    GvsBoundBox   operator*  ( const GvsBoundBox &otherBox ) const;   // intersection
    GvsBoundBox&  operator*= ( const GvsBoundBox &otherBox );

    bool          getTentryTexit ( const m4d::vec3 &p0, const m4d::vec3 &p1, double tp0, double tp1,
                                   double &time_Entry, double &time_Exit,
                                   short &entryFace, short &exitFace ) const;

    virtual   void  scale      ( const m4d::vec3&   scaleVec );
    virtual   void  translate  ( const m4d::vec3&   transVec );
    virtual   void  rotate     ( const m4d::vec3&   rotAxis, double rotAngle );
    virtual   void  transform  ( const m4d::Matrix<double,3,4>& mat      );

    int       getChart ( ) const;
    void      setChart ( const int chart );

    void      Print ( FILE* fptr = stderr );

  protected:
    m4d::vec3  boxLowBounds;
    m4d::vec3  boxUppBounds;
};


inline
m4d::vec3 GvsBoundBox :: lowBounds ( ) const
{
  return boxLowBounds;
}

inline
m4d::vec3 GvsBoundBox :: uppBounds ( ) const
{
  return boxUppBounds;
}

#endif
