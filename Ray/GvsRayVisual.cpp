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

#include "Ray/GvsRayVisual.h"


double GvsRayVisual :: rayContribThreshold = 0.01;
ushort GvsRayVisual :: rayMaxTreeDepth = 0;


static int getRayTreeDepth ( int rayIndex )
{
    if (rayIndex == 0)
        return 0;

    if (rayIndex < 0 )
        return -rayIndex;

    int rayTreeDepth = 1;

    while ( (rayIndex >>= 1) > 0 )
        rayTreeDepth += 1;

    return rayTreeDepth;
}

GvsRayVisual :: GvsRayVisual ( ) : GvsRayClosestIS () {
}

GvsRayVisual :: GvsRayVisual ( GvsRayGen* gen )
    : GvsRayClosestIS ( gen )
{
}

GvsRayVisual :: GvsRayVisual ( const m4d::vec4 &orig, const m4d::vec4 &dir, GvsRayGen* gen, int index )
    : GvsRayClosestIS ( orig, dir, gen )
{
    rayIndex = index;
    rayTreeDepth = getRayTreeDepth ( index );
    rayMaxContrib = 1.0;
}

GvsRayVisual :: GvsRayVisual ( const m4d::vec4 &orig, const m4d::vec4 &dir, GvsRayGen* gen,
                               double minSearchDist, double maxSearchDist, int index )
    : GvsRayClosestIS ( orig, dir, gen, minSearchDist, maxSearchDist )
{
    rayIndex = index;
    rayTreeDepth = getRayTreeDepth ( index );
    rayMaxContrib = 1.0;
}

GvsRayVisual :: GvsRayVisual ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad,
                               GvsRayGen* gen, int index )
    : GvsRayClosestIS ( orig, dir, tetrad, gen )
{
    rayIndex = index;
    rayTreeDepth = getRayTreeDepth ( index );
    rayMaxContrib = 1.0;
}

GvsRayVisual :: GvsRayVisual ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad,
                               GvsRayGen* gen, double minSearchDist, double maxSearchDist, int index )
    : GvsRayClosestIS ( orig, dir, tetrad, gen, minSearchDist, maxSearchDist )
{
    rayIndex = index;
    rayTreeDepth = getRayTreeDepth ( index );
    rayMaxContrib = 1.0;
}

GvsRayVisual :: GvsRayVisual ( m4d::vec4* geodPoints, m4d::vec4* geodTangents, int numPoints, int index )
    : GvsRayClosestIS ( )
{
    setPoints( geodPoints );
    setDirs( geodTangents );
    setNumPoints( numPoints );

    rayIndex = index;
    rayTreeDepth = getRayTreeDepth ( index );
    rayMaxContrib = 1.0;
    //cout << minSearchDist() << " " << maxSearchDist() << std::endl;
}

GvsRayVisual::~GvsRayVisual() {
}

m4d::vec3 GvsRayVisual :: getReflectedRayDir() {
    m4d::vec3 rayDirection = raySurfIntersec.getLocalDirection();
    //rayDirection.print(cerr);
    assert ( rayDirection!=m4d::vec3() );

    m4d::vec3 reflecRayDir = rayDirection - 2.0* (raySurfIntersec.normal() | rayDirection) * raySurfIntersec.normal();
    reflecRayDir.normalize();

    return reflecRayDir;
}


GvsRayVisual GvsRayVisual :: getReflectedRay() {
    m4d::vec3 reflecRayDir = getReflectedRayDir();

    GvsRayVisual reflecRay(raySurfIntersec.point(),reflecRayDir.get3As4(),getRayGen());

    reflecRay.rayTreeDepth = rayTreeDepth + 1;
    reflecRay.rayIndex = rayIndex + rayIndex + 1;

    return reflecRay;
}


GvsRayVisual GvsRayVisual :: getTransmittedRay ( double  ) {
    std::cerr << "GvsRayVisual :: getTransmittedRay() ... not implemented yet!\n";
    return GvsRayVisual();
}
