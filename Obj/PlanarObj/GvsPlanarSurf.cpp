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

#include "Obj/PlanarObj/GvsPlanarSurf.h"

#include <Ray/GvsRay.h>
#include <Ray/GvsRayAllIS.h>


GvsPlanarSurf :: GvsPlanarSurf ( GvsSurfaceShader* shader )
    :  GvsSurface ( shader )
{
    planeNormal = m4d::vec3(0,0,1);
    planeDist = 0.0;    
    mID = ++mObjCounter;
}


GvsPlanarSurf :: GvsPlanarSurf ( const GvsPlanarSurf& surf )
    :  GvsSurface ( surf.surfShader )
{
    planeNormal = m4d::vec3(0,0,1);
    planeDist = 0.0;
    mID = ++mObjCounter;
}


GvsPlanarSurf :: ~GvsPlanarSurf() {
}

m4d::vec3 GvsPlanarSurf :: normal() const {
    return planeNormal;
}

bool GvsPlanarSurf :: testIntersection ( GvsRay& ray ) {    
    assert( mMetric != NULL );

    int chart0,chart1;

    int maxSeg   = ray.getNumPoints()-2;
    int startSeg = int(0);
    int endSeg   = maxSeg;

    m4d::enum_coordinate_type coords = mMetric->getCoordType();

    for( int seg=startSeg; seg <= endSeg; seg++)  {
        m4d::vec4 p0 = ray.getPoint(seg);
        m4d::vec4 p1 = ray.getPoint(seg+1);

        m4d::vec4 p0trans4D = p0;
        m4d::vec4 p1trans4D = p1;

        chart0 = chart1 = 0;
        if (coords != m4d::enum_coordinate_cartesian)  {
            chart0 = mMetric->transToPseudoCart( p0, p0trans4D );
            chart1 = mMetric->transToPseudoCart( p1, p1trans4D );
        }

        if (chart0!=mChart && chart1!=mChart) {
            continue;
        }

        if (haveMotion) {
            m4d::vec4 p0motTrans4D = p0trans4D;
            m4d::vec4 p1motTrans4D = p1trans4D;
            stMotion->getTransformedPolygon(seg,p0motTrans4D,p1motTrans4D, p0trans4D, p1trans4D);
        }
        m4d::vec4 vtrans4D  = p1trans4D - p0trans4D;

        m4d::vec3 p0trans = p0trans4D.getAsV3D();
        m4d::vec3 p1trans = p1trans4D.getAsV3D();
        m4d::vec3 vtrans = p1trans - p0trans;

        double tp0 = p0.x(0);
        double tp1 = p1.x(0);

        double tHit,alpha;
        m4d::vec3 rayIntersecPt;
        if (rayIntersect(p0trans,p1trans,tp0,tp1,alpha,tHit,rayIntersecPt)) {
            if (GvsRay::isIn(seg,alpha,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,alpha))) {
                if (isValidHit(rayIntersecPt)) {
                    //std::cerr << "intersec " << alpha << " " << tHit << std::endl;

                    GvsSurfIntersec surfIntersec;
                    surfIntersec.setDist ( GvsRay::calcRayDist(seg,alpha) );
                    surfIntersec.setSurface ( this );

                    // global intersection point and direction in proper metric coordinates
                    m4d::vec4 point = p0trans4D + alpha * vtrans4D;
                    m4d::vec4 dir   = vtrans4D;
                    m4d::enum_coordinate_type cType = mMetric->getCoordType();
                    m4d::TransCoordinates::coordTransf(m4d::enum_coordinate_cartesian,point,dir,cType,point,dir);
                    surfIntersec.setPoint( point );
                    surfIntersec.setDirection( dir );

                    // local intersection point in standard object system
                    surfIntersec.setLocalPoint( p0trans + alpha * vtrans );
                    surfIntersec.setLocalDirection( vtrans );

                    surfIntersec.setRaySegNumber(seg);
                    return ray.store( surfIntersec );
                }
            }
        }

        // TODO

    }
    return false;
}


bool GvsPlanarSurf::isValidHit ( m4d::vec3 ) {
    return true;
}

bool GvsPlanarSurf :: allIntersections ( GvsRayAllIS& ray ) {
    return testIntersection( ray );
}


bool GvsPlanarSurf :: testLocalIntersection ( GvsRay&     ,
                                              const int ,
                                              const int ,
                                              const m4d::vec4 ,
                                              const m4d::vec4  )
{
    //pure virtual
    return false;
}


bool
GvsPlanarSurf :: testLocalIntersection ( GvsRay&     ,
                                         const int ,
                                         const int ,
                                         const m4d::vec4 ,
                                         const m4d::vec4 ,
                                         GvsSurface*  )
{
    //pure virtual
    return false;
}


bool
GvsPlanarSurf :: allLocalIntersections ( GvsRayAllIS& ray,
                                         const int rayPartIndex,
                                         const int seg,
                                         const m4d::vec4 p0,
                                         const m4d::vec4 p1 )
{
    return testLocalIntersection( ray, rayPartIndex,seg, p0, p1 );
}


GvsBoundBox
GvsPlanarSurf :: boundingBox ( void ) const
{
    return planarSurfBoundBox;
}



bool GvsPlanarSurf :: rayIntersect (const m4d::vec3& p0, const m4d::vec3& p1,
                                    double tp0, double tp1, double &alpha,
                                    double &thit,
                                    m4d::vec3& rayIntersecPnt ) const
{
    double delta_t = tp1-tp0;

    m4d::vec3 rayOrig3D = p0;
    m4d::vec3 rayDir3D  = p1-p0;

    double cosRayPlane = planeNormal | rayDir3D;
    if ( fabs(cosRayPlane) > GVS_EPS ) {
        alpha = -((planeNormal|rayOrig3D) - planeDist)/cosRayPlane;
        //fprintf(stderr," %f %f\n",alpha,planeDist);
        thit = tp0 + delta_t * alpha;
        rayIntersecPnt = rayOrig3D + alpha*rayDir3D;
        return true;
    }
    return false;
}

