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
#ifndef GVS_RAY_H
#define GVS_RAY_H

#include "GvsGlobalDefs.h"
#include "Obj/STMotion/GvsLocalTetrad.h"
//#include <Shader/GvsShader.h>
//#include <Shader/Surface/GvsSurfaceShader.h>

#include <metric/m4dMetric.h>

#include <cassert>
#include <vector>
#include <iostream>
#include <limits.h>


enum GvsRayType {
    polRay, polRayOneIS, polRayClosestIS, polRayVisual,
    polRayAnyIS, polRayAllIS
};

class GvsRayGen;
class GvsSurfIntersec;


class GvsRay
{
public:
    GvsRay ( );
    GvsRay ( GvsRayGen* gen );
    GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, GvsRayGen* gen );
    GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, GvsRayGen* gen,
             double minSearchDist, double maxSearchDist );

    GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad, GvsRayGen* gen );
    GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad, GvsRayGen* gen,
             double minSearchDist, double maxSearchDist );

    virtual ~GvsRay();

    
    virtual bool   recalc ( const m4d::vec4 &orig, const m4d::vec4 &dir );
    virtual bool   recalc ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad* tetrad );
    virtual bool   recalcJacobi ( const m4d::vec4 &orig, const m4d::vec4 &dir,
                                  const m4d::vec3 &locRayDir, const GvsLocalTetrad* tetrad );

    void           setSearchInterval ( double minDist, double maxDist );

    GvsRayGen*     getRayGen    () const;
    m4d::vec4*     points       ();
    m4d::vec4*     tangents     ();
    int            getNumPoints () const;

    void         setPoints      ( m4d::vec4* pts  );
    void         setDirs        ( m4d::vec4* dirs );
    void         setNumPoints   ( int noPts );

    m4d::vec4      getPoint       ( int index ) const;
    m4d::vec4      getTangente    ( int index ) const;
    m4d::vec5      getJacobi      ( int index ) const;
    GvsLocalTetrad getTetrad      ( int index ) const;


    ulong          getID          ( ) const;

    double         minSearchDist  ( ) const;
    double         maxSearchDist  ( ) const;

    virtual bool   intersecFound  ( ) const = 0;
    virtual bool   store          ( const GvsSurfIntersec &surfIntersec ) = 0;  // pure virtual
    virtual bool   isValidSurfIntersec ( double dist ) const;

    virtual void   Print ( FILE* fptr = stderr );



    static ulong   getNonRayID       ( void );
    static ulong   getNumRaysTraced  ( void );

    void           timeShiftRay      ( double timeDelta );

    static double  calcRayDist  ( const int seg, const double alpha );
    static bool    isIn         ( const int seg, const double alpha , const double endDist);
    static bool    isValidAlpha ( const double alpha );
    static void    splitRayDist ( const double distance, int &seg, double &alpha );


protected:
    void  setMinSearchDist ( double minDist );
    void  setMaxSearchDist ( double maxDist );

    void deleteAll();

private:
    static ulong  getNextRayID();


private:
    ulong               rayID;
    GvsRayGen*          rayGen;
    m4d::vec4*          rayPoints;
    m4d::vec4*          rayDirs;
    GvsLocalTetrad*     rayTetrad;    
    double*             rayLambda;
    m4d::vec4*          raySachs1;
    m4d::vec4*          raySachs2;
    m4d::vec5*          rayJacobi;
    m4d::vec5           rayMaxJacobi;

    int                 rayNumPoints;
    double              rayMinSearchDist;
    double              rayMaxSearchDist;

    bool                rayHasTetrad;
};

//----------------------------------------------------------------------------
//            i n l i n e
//----------------------------------------------------------------------------
inline ulong GvsRay :: getNonRayID () {
    return 0UL;
}

inline ulong GvsRay :: getNumRaysTraced () {
  return getNextRayID() - 1;
}

#endif
