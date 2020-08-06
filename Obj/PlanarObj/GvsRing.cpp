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

#include "Obj/PlanarObj/GvsRing.h"

#include <Ray/GvsRay.h>
#include <Ray/GvsRayAllIS.h>
#include <Obj/STMotion/GvsStMotion.h>

#include "math/TransfMat.h"


GvsRing::GvsRing()
    :  GvsPlanarSurf( NULL )
{
    ringCenter = m4d::vec3(0.0,0.0,0.0);
    ringNormal = m4d::vec3(0.0,0.0,1.0);
    planeDist = 0.0;

    ringOuterRadius = 2.0;
    ringInnerRadius = 1.0;

    mMetric = NULL;
    stMotion = NULL;
    mObjType = inCoords;
    calcBoundBox();
}


GvsRing :: GvsRing(const m4d::vec3& center, const m4d::vec3& normal,
                   const double rOuter, const double rInner,
                   GvsSurfaceShader* shader, m4d::Metric* metric,
                   GvsObjType objType )
    :  GvsPlanarSurf ( shader )
{
    ringCenter = center;
    ringNormal = normal.getNormalized();

    // TODO: planeDist = ?

    ringOuterRadius = rOuter;
    ringInnerRadius = rInner;

    mMetric = metric;
    stMotion = NULL;
    mObjType = objType;
    calcBoundBox();
}

GvsRing::~GvsRing(){
}


GvsPlanarSurf* GvsRing :: getClone( void ) const {
    GvsRing* clone = new GvsRing( *this );
    assert( clone != NULL );
    return clone;
}

bool GvsRing::isValidHit ( m4d::vec3 rp ) {
    double dist = (rp-ringCenter).getNorm();
    if (dist>=ringInnerRadius && dist<=ringOuterRadius) {
        return true;
    }
    return false;
}


void GvsRing::calcNormal( GvsSurfIntersec & intersec ) const {
    intersec.setNormal( ringNormal );
}


void GvsRing::calcTexUVParam( GvsSurfIntersec & intersec ) const {
    double u = intersec.localPoint().x(0);
    double v = intersec.localPoint().x(1);

    double r = (sqrt(u*u+v*v)-ringInnerRadius)/(ringOuterRadius-ringInnerRadius);
    double phi = atan2(v,u)/(2.0*M_PI) + 0.5;
    //fprintf(stderr," %f %f  %f %f\n",u,v,r,phi);
    intersec.setTexUVParam(phi,r);
}


void GvsRing::Print ( FILE* fptr ) {
    fprintf(fptr,"Ring {\n");
    fprintf(fptr,"\tcenter: ");ringCenter.printS(fptr);
    fprintf(fptr,"\tnormal: ");ringNormal.printS(fptr);
    fprintf(fptr,"\tRout:   %f\n",ringOuterRadius);
    fprintf(fptr,"\tRin:    %f\n",ringInnerRadius);
    fprintf(fptr,"}\n");
}


void GvsRing::calcBoundBox() {
    m4d::vec3 z = m4d::vec3(0.0,0.0,1.0);
    m4d::vec3 v1 = ringNormal^z;
    if (v1.isZero()) {
        z = m4d::vec3(1.0,0.0,0.0);
        v1 = ringNormal^z;
    }
    v1 = v1.getNormalized();

    m4d::vec3 v2 = ringNormal^v1;
    planarSurfBoundBox = GvsBoundBox( ringCenter+ringOuterRadius*v1, ringCenter+ringOuterRadius*v2 );
}


