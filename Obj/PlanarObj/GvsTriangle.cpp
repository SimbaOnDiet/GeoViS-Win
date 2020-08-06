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

#include "Obj/PlanarObj/GvsTriangle.h"

#include <Ray/GvsRay.h>
#include <Ray/GvsRayAllIS.h>
#include <Obj/STMotion/GvsStMotion.h>

#include <cassert>

#include "math/TransfMat.h"

double area( const m4d::vec2 &v1, const m4d::vec2 &v2 ) {
    return v1.x(0) * v2.x(1) - v1.x(1) * v2.x(0);
}

double area( const m4d::vec3& v1, const m4d::vec3& v2 ) {
    return (v1^v2).getNorm();
}

//------------------------------------------------------------------------
//    constructor  	G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: GvsTriangle 	( void )
    :  GvsPlanarSurf( NULL )
{
    // STMotion ??
    trianVertex[0] = m4d::vec3( 0.0, 0.0, 0.0 );
    trianVertex[1] = m4d::vec3( 1.0, 0.0, 0.0 );
    trianVertex[2] = m4d::vec3( 0.0, 1.0, 0.0 );

    trianUV[0] = m4d::vec2(0.0,0.0);
    trianUV[1] = m4d::vec2(1.0,0.0);
    trianUV[2] = m4d::vec2(0.0,1.0);

    planeNormal = m4d::vec3(0.0,0.0,1.0);

    mMetric = NULL;
    stMotion = NULL;
    mObjType = inCoords;
    calcBoundBox();
}


//------------------------------------------------------------------------
//    constructor  	G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: GvsTriangle 	( const m4d::vec3& p0, const m4d::vec3& p1,
                              const m4d::vec3& p2, GvsSurfaceShader* shader,
                              m4d::Metric* metric, GvsObjType objType)
    :  GvsPlanarSurf ( shader )
{
    trianVertex[0] = p0;
    trianVertex[1] = p1;
    trianVertex[2] = p2;

    planeNormal = calcTriNormal();

    trianUV[0] = m4d::vec2(0.0,0.0);
    trianUV[1] = m4d::vec2(1.0,0.0);
    trianUV[2] = m4d::vec2(0.0,1.0);

    mMetric = metric;
    stMotion = NULL;
    mObjType = objType;
    calcBoundBox();
}


//------------------------------------------------------------------------
//    constructor  	G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: GvsTriangle 	( const m4d::vec3& p0, const m4d::vec3& p1,
                              const m4d::vec3& p2, GvsSurfaceShader* shader,
                              m4d::Metric* metric, GvsStMotion* motion,
                              GvsObjType objType )
    :  GvsPlanarSurf ( shader )
{
    trianVertex[0] = p0;
    trianVertex[1] = p1;
    trianVertex[2] = p2;

    planeNormal = calcTriNormal();

    trianUV[0] = m4d::vec2(0.0,0.0);
    trianUV[1] = m4d::vec2(1.0,0.0);
    trianUV[2] = m4d::vec2(0.0,1.0);

    mMetric = metric;
    stMotion = motion;
    mObjType = objType;
    calcBoundBox();
}

//------------------------------------------------------------------------
//    constructor  	G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: GvsTriangle 	( const m4d::vec3& p0,  const m4d::vec3& p1,
                              const m4d::vec3& p2,  const m4d::vec2& uv0,
                              const m4d::vec2& uv1, const m4d::vec2& uv2,
                              GvsSurfaceShader* shader,
                              m4d::Metric* metric, GvsObjType objType)
    :  GvsPlanarSurf ( shader )
{
    trianVertex[0] = p0;
    trianVertex[1] = p1;
    trianVertex[2] = p2;

    planeNormal = calcTriNormal();

    trianUV[0] = uv0;
    trianUV[1] = uv1;
    trianUV[2] = uv2;

    mMetric = metric;
    stMotion = NULL;
    mObjType = objType;
    calcBoundBox();
}


//------------------------------------------------------------------------
//    constructor  	G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: GvsTriangle 	( const m4d::vec3& p0,  const m4d::vec3& p1,
                              const m4d::vec3& p2,  const m4d::vec2& uv0,
                              const m4d::vec2& uv1, const m4d::vec2& uv2,
                              GvsSurfaceShader* shader,
                              m4d::Metric* metric, GvsStMotion* motion,
                              GvsObjType objType )
    :  GvsPlanarSurf ( shader )
{
    trianVertex[0] = p0;
    trianVertex[1] = p1;
    trianVertex[2] = p2;

    planeNormal = calcTriNormal();

    trianUV[0] = uv0;
    trianUV[1] = uv1;
    trianUV[2] = uv2;

    mMetric = metric;
    stMotion = motion;
    mObjType = objType;
}



//------------------------------------------------------------------------
//    copyconstructor  	G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: GvsTriangle 	( const GvsTriangle& triangle )
    :  GvsPlanarSurf ( triangle )
{
    trianVertex[0] = triangle.trianVertex[0];
    trianVertex[1] = triangle.trianVertex[1];
    trianVertex[2] = triangle.trianVertex[2];

    planeNormal = triangle.planeNormal;

    trianUV[0] = triangle.trianUV[0];
    trianUV[1] = triangle.trianUV[1];
    trianUV[2] = triangle.trianUV[2];

    mMetric  = triangle.mMetric;
    stMotion = triangle.stMotion;
    mObjType = triangle.mObjType;
    calcBoundBox();
}


//------------------------------------------------------------------------
//    destructor   	~ G v s T r i a n g l e
//------------------------------------------------------------------------
GvsTriangle :: ~GvsTriangle 	( void )
{

}


//------------------------------------------------------------------------
//     method		g e t C l o n e
//------------------------------------------------------------------------
GvsPlanarSurf*
GvsTriangle :: getClone		( void ) const
{
    GvsTriangle* clone = new GvsTriangle( *this );
    assert( clone != NULL );
    return clone;
}

bool GvsTriangle::isValidHit ( m4d::vec3 rp ) {
    m4d::vec3 vq1 = trianVertex[0]-rp;
    m4d::vec3 vq2 = trianVertex[1]-rp;
    m4d::vec3 vq3 = trianVertex[2]-rp;

    double Area12 = area(vq1,vq2);
    double Area13 = area(vq1,vq3);
    double Area23 = area(vq2,vq3);
    //fprintf(stderr," is: %f %f\n",Area12+Area13+Area23,trianArea);
    if (Area12+Area13+Area23 <= 2.0*trianArea+GVS_EPS) {
        return true;
    }
    return false;
}


m4d::vec3 GvsTriangle::calcTriNormal() {
    m4d::vec3 p1 = trianVertex[0];
    m4d::vec3 p2 = trianVertex[1];
    m4d::vec3 p3 = trianVertex[2];
    m4d::vec3 n = (p2-p1)^(p3-p1);
    trianArea = 0.5*n.getNorm();
    n = n.getNormalized();
    planeDist = p1|n;
    return n;
}

//------------------------------------------------------------------------
//     method		c a l c N o r m a l
//------------------------------------------------------------------------
void GvsTriangle :: calcNormal( GvsSurfIntersec & intersec ) const {
    intersec.setNormal( planeNormal );
}


//------------------------------------------------------------------------
//     method		c a l c T e x U V P a r a m
//------------------------------------------------------------------------
void
GvsTriangle :: calcTexUVParam	( GvsSurfIntersec & intersec ) const
{
    /* F?r die Default-Werte von trianUV ergibt sich:
     intersec.setTexUVParam( v, w )
  */

    double u = intersec.surfSTParam().x(0);
    double v = intersec.surfSTParam().x(1);
    double w = 1.0 - u - v;

    intersec.setTexUVParam(
                u * trianUV[0].x(0) + v * trianUV[1].x(0) + w * trianUV[2].x(0),
            u * trianUV[0].x(1) + v * trianUV[1].x(1) + w * trianUV[2].x(1)  );
}


//------------------------------------------------------------------------
//     method		s c a l e
//------------------------------------------------------------------------
void
GvsTriangle :: scale ( const m4d::vec3&   scaleVec )
{
    for ( int i = 0; i < 3; i++ )
    {
        trianVertex[i] *= scaleVec;
    }

    planeNormal = calcTriNormal();
    calcBoundBox();
}


//------------------------------------------------------------------------
//     method		t r a n s l a t e
//------------------------------------------------------------------------
void
GvsTriangle :: translate( const m4d::vec3&   transVec )
{
    for ( int i = 0; i < 3; i++ ) {
        trianVertex[i] += transVec;
    }
    calcBoundBox();
}


//------------------------------------------------------------------------
//     method		r o t a t e
//------------------------------------------------------------------------
void
GvsTriangle :: rotate ( const m4d::vec3&   rotAxis,
                        const double rotAngle )
{
    m4d::RotateMat3D rotMat = m4d::RotateMat3D( rotAxis, rotAngle );
    for ( int i = 0; i < 3; i++ ) {
        trianVertex[i] = rotMat * trianVertex[i];
    }
    planeNormal = rotMat * planeNormal;
}


void GvsTriangle :: transform( const m4d::Matrix<double,3,4>& mat ) {
    for ( int i = 0; i < 3; i++ ) {
        //trianVertex[i] = mat * trianVertex[i];
    }

    //planeNormal = rotMat * planeNormal;
    calcBoundBox();
}


//------------------------------------------------------------------------
//     method  p r i n t
//------------------------------------------------------------------------
void GvsTriangle::Print ( FILE* fptr ) {
    fprintf(stderr,"Triangle {\n");
    fprintf(stderr,"\tv1: ");trianVertex[0].printS(fptr);
    fprintf(stderr,"\tv2: ");trianVertex[1].printS(fptr);
    fprintf(stderr,"\tv3: ");trianVertex[2].printS(fptr);
    fprintf(stderr,"}\n");
}


//------------------------------------------------------------------------
//     method 		c a l c P r o j E d g e s
//------------------------------------------------------------------------
void
GvsTriangle :: calcProjEdges	( int          projPlane,
                                  const m4d::vec3&   intersecPt,
                                  m4d::vec2&         projEdge20,
                                  m4d::vec2&         projEdge21,
                                  m4d::vec2&         projEdge2p ) const
{
    switch ( projPlane )
    {
        case 0:
            projEdge20.x(0) = trianVertex[0].x(1) - trianVertex[2].x(1);
            projEdge20.x(1) = trianVertex[0].x(2) - trianVertex[2].x(2);
            projEdge21.x(0) = trianVertex[1].x(1) - trianVertex[2].x(1);
            projEdge21.x(1) = trianVertex[1].x(2) - trianVertex[2].x(2);
            projEdge2p.x(0) = intersecPt.x(1) - trianVertex[2].x(1);
            projEdge2p.x(1) = intersecPt.x(2) - trianVertex[2].x(2);
            break;
        case 1:
            projEdge20.x(0) = trianVertex[0].x(0) - trianVertex[2].x(0);
            projEdge20.x(1) = trianVertex[0].x(2) - trianVertex[2].x(2);
            projEdge21.x(0) = trianVertex[1].x(0) - trianVertex[2].x(0);
            projEdge21.x(1) = trianVertex[1].x(2) - trianVertex[2].x(2);
            projEdge2p.x(0) = intersecPt.x(0) - trianVertex[2].x(0);
            projEdge2p.x(1) = intersecPt.x(2) - trianVertex[2].x(2);
            break;
        case 2:
            projEdge20.x(0) = trianVertex[0].x(0) - trianVertex[2].x(0);
            projEdge20.x(1) = trianVertex[0].x(1) - trianVertex[2].x(1);
            projEdge21.x(0) = trianVertex[1].x(0) - trianVertex[2].x(0);
            projEdge21.x(1) = trianVertex[1].x(1) - trianVertex[2].x(1);
            projEdge2p.x(0) = intersecPt.x(0) - trianVertex[2].x(0);
            projEdge2p.x(1) = intersecPt.x(1) - trianVertex[2].x(1);
            break;
    }
}


//------------------------------------------------------------------------
//     method		c a l c B o u n d B o x
//------------------------------------------------------------------------
void
GvsTriangle :: calcBoundBox	( void )
{
    planarSurfBoundBox = GvsBoundBox ( trianVertex[0], trianVertex[1] );
    planarSurfBoundBox.extendBoxToContain ( trianVertex[2] );
}


