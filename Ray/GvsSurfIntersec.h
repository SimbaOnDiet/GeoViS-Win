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
#ifndef GVS_SURF_INTERSEC_H
#define GVS_SURF_INTERSEC_H

#include "GvsGlobalDefs.h"

#include "Obj/PlanarObj/GvsPlanarSurf.h"
#include "Obj/STMotion/GvsLocalTetrad.h"
#include "Shader/GvsShader.h"


class GvsSurface;

class GvsSurfIntersec
{
public:
    enum InsecSubType { Default, Entering, Exiting };

public:
    GvsSurfIntersec ();
    GvsSurfIntersec ( const GvsSurfIntersec& surfIntersec );
    virtual ~GvsSurfIntersec ();

    GvsSurfIntersec& operator = ( const GvsSurfIntersec& surfIntersec );

    virtual double       dist   () const;
    virtual GvsShader*   shader () const;
    virtual GvsSceneObj* object () const;
    virtual void         reset  ();


    // ---------------------------------------------------

    m4d::vec4     point                () const;
    m4d::vec4     direction            () const;
    InsecSubType  subType              () const;
    GvsSurface*   surface              () const;
    bool          surfIsSelfDescribing () const;

    m4d::vec3     normal      ();
    m4d::vec3     normNormal  ();
    m4d::vec3     derivS      ();
    m4d::vec3     derivT      ();
    m4d::vec2     texUVParam  ();

    m4d::vec2     surfSTParam  () const;
    m4d::vec3     localPoint   () const;

    bool      surfSTParamAreValid ( void                ) const;
    bool      localPointIsValid   ( void                ) const;

    // ---------------------------------------------------
    void      setDist             ( double              );
    void      setPoint            ( const m4d::vec4 &point );
    void      setDirection        ( const m4d::vec4 &dir   );
    void      setSubType          ( const InsecSubType& );
    void      setSurface          ( GvsSurface*         );

    void      setSurfIsSelfDescribing ( bool            );

    void      setInstSurfIntersec ( GvsSurfIntersec*    );
    void      setTex3dPoint       ( const m4d::vec3& );
    void      setNormal           ( const m4d::vec3& );
    void      setDerivS           ( const m4d::vec3& );
    void      setDerivT           ( const m4d::vec3& );
    void      setTexUVParam       ( const m4d::vec2& );
    void      setTexUVParam       ( double, double );

    void      setSurfSTParam      ( const m4d::vec2&  );
    void      setSurfSTParam      ( double, double );
    void      setLocalPoint       ( const m4d::vec3&  );

    void      invertNormal        ();

    void      setRaySegNumber     ( int segNumber );
    int       getRaySegNumber     ();

    void            setLocalIntersec ( bool locIntersec   );
    bool            isLocalIntersec  () const;

    void            setLocalTetrad   ( GvsLocalTetrad *lt );
    GvsLocalTetrad* getLocalTetrad   () const;


    void       setLocalDirection   ( const m4d::vec3& localDir );
    m4d::vec3  getLocalDirection   () const;
    
    void       setLocalTime        ( double localTime );
    double     getLocalTime        () const;

    void print ( FILE* fptr = stderr );


    bool operator == ( double rayLength    ) const;
    bool operator != ( double rayLength    ) const;
    bool operator <  ( double rayLength    ) const;
    bool operator <= ( double rayLength    ) const;
    bool operator >  ( double rayLength    ) const;
    bool operator >= ( double rayLength    ) const;

    bool operator == ( const GvsSurfIntersec & ) const;
    bool operator != ( const GvsSurfIntersec & ) const;
    bool operator <  ( const GvsSurfIntersec & ) const;
    bool operator <= ( const GvsSurfIntersec & ) const;
    bool operator >  ( const GvsSurfIntersec & ) const;
    bool operator >= ( const GvsSurfIntersec & ) const;

private:
    void assign( const GvsSurfIntersec& );

private:
    // Intersection parameters that have to be determined every time
    double         insecDist;
    m4d::vec4      insecPoint;
    m4d::vec4      insecDirection;
    InsecSubType   insecSubType;
    GvsSurface*    insecSurface;
    bool           insecSurfIsSelfDescribing;
    int            insecRaySegNumber;

    // Intersection parameters determined by callback.
    m4d::vec3      insecNormal;
    m4d::vec3      insecNormNormal;
    bool           insecNormalIsValid;
    m4d::vec3      insecDerivS;
    bool           insecDerivSIsValid;
    m4d::vec3      insecDerivT;
    bool           insecDerivTIsValid;
    m4d::vec2      insecTexUVParam;
    bool           insecTexUVParamAreValid;

    // Optional intersection parameters.
    m4d::vec3      insecLocalPoint;
    bool           insecLocalPointIsValid;
    m4d::vec2      insecSurfSTParam;
    bool           insecSurfSTParamAreValid;
    m4d::vec3      insecLocalDir;         // Direction of incoming light ray (local)
    double         insecLocalTime;        // proper time

    bool             insecLocalInsec;
    GvsLocalTetrad*  insecLocalTetrad;

public:
    int   partIndex;    // part is usually a face (box)
};


#endif
