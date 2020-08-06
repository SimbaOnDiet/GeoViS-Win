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

#include "Obj/SolidObj/GvsSolEllipsoid.h"
#include "Obj/SolidObj/GvsSolObjSpanList.h"
#include "Ray/GvsRay.h"
#include "Ray/GvsSurfIntersec.h"
#include "Shader/Surface/GvsSurfaceShader.h"

#include <metric/m4dMetric.h>
#include "math/TransfMat.h"


GvsSolEllipsoid :: GvsSolEllipsoid( const m4d::vec3&  center,
                                    const m4d::vec3&  halfAxisLength,
                                    GvsSurfaceShader* shader,
                                    m4d::Metric*      metric,
                                    GvsObjType        objType )
    : GvsSolConvexPrim(shader)
{
    mObjType = objType;
    mMetric = metric;
    setGeometry( center, halfAxisLength );
}




GvsSolEllipsoid :: GvsSolEllipsoid( const m4d::vec3&  center,
                                    const m4d::vec3&  halfAxisLength,
                                    GvsSurfaceShader* shader,
                                    m4d::Metric*      metric,
                                    GvsStMotion*      motion,
                                    GvsObjType        objType )
    : GvsSolConvexPrim( shader )
    
{
    mObjType = objType;
    mMetric = metric;
    stMotion = motion;
    setGeometry( center, halfAxisLength );
}


void GvsSolEllipsoid :: setGeometry(const m4d::vec3& center, const m4d::vec3& halfAxisLength) {
    setCenter( center );
    setHalfAxisLength( halfAxisLength );
    
    volTransfMat    = m4d::TranslateMat3D(center) * m4d::ScaleMat3D(halfAxisLength);
    volInvTransfMat = m4d::ScaleMat3D(1.0/ halfAxisLength.x(0),
                                      1.0/ halfAxisLength.x(1),
                                      1.0/ halfAxisLength.x(2)) *
            m4d::TranslateMat3D(-center);
    calcBoundBox();
}


void GvsSolEllipsoid :: setCenter( m4d::vec3 center ) {
    mCenter = center;
}

void GvsSolEllipsoid :: setHalfAxisLength( m4d::vec3 halfAxisLength ) {
    mHalfAxisLength = halfAxisLength;
}


bool GvsSolEllipsoid::getTentryTexit( const m4d::vec3& p0, const m4d::vec3& p1,
                                      double tp0, double tp1,
                                      double& time_Entry, double& time_Exit,
                                      short&, short&) const
{
    double delta_tp = tp1 - tp0;

    double mx = ( p1.x(0) - p0.x(0) ) / delta_tp;
    double my = ( p1.x(1) - p0.x(1) ) / delta_tp;
    double mz = ( p1.x(2) - p0.x(2) ) / delta_tp;

    double ax = p0.x(0) - mx*tp0;
    double ay = p0.x(1) - my*tp0;
    double az = p0.x(2) - mz*tp0;

    double a = mx*mx + my*my + mz*mz;
    double b = 2.0*( ax*mx + ay*my + az*mz);
    double c = ax*ax + ay*ay + az*az - 1.0;

    double d = b*b - ( 4.0*a*c);

    if ( d < 0.0) {
        return false;
    }
    d = sqrt(d);

    time_Entry = ( - b - d ) / (2.0 * a);
    time_Exit  = ( - b + d ) / (2.0 * a);

    if ( tp0 < tp1 ) {
        if (time_Entry > time_Exit) {
            double t   = time_Entry;
            time_Entry = time_Exit;
            time_Exit  = t;
        }
    }
    else {
        if (time_Entry < time_Exit) {
            double t   = time_Entry;
            time_Entry = time_Exit;
            time_Exit  = t;
        }
    }
    return true;
}


m4d::vec3 GvsSolEllipsoid :: getCenter() const {
    return mCenter;
}

m4d::vec3 GvsSolEllipsoid :: getHalfAxisLength() const {
    return mHalfAxisLength;
}


bool GvsSolEllipsoid::getRaySpanList( GvsRay& ray, GvsSolObjSpanList& spanList ) {
    double time_Entry, time_Exit;
    double tEntry, tExit;
    int chart0, chart1;

    bool entryFound = false;
    bool exitFound  = false;

    GvsSurfIntersec insecEntry, insecExit;

    m4d::Metric* stMetric = mMetric;
    m4d::enum_coordinate_type coords = stMetric->getCoordType();

    int maxSeg   = ray.getNumPoints()-2;
    int startSeg = int(0);
    int endSeg   = maxSeg;


    for( int seg = startSeg; seg < endSeg; seg++ ) {
        m4d::vec4 p0 = ray.getPoint(seg);
        m4d::vec4 p1 = ray.getPoint(seg+1);

        m4d::vec4 p0trans4D = p0;
        m4d::vec4 p1trans4D = p1;

        chart0 = chart1 = 0;
        if (coords != m4d::enum_coordinate_cartesian) {
            chart0 = stMetric->transToPseudoCart( p0, p0trans4D );
            chart1 = stMetric->transToPseudoCart( p1, p1trans4D );
        }

        if (chart0!=mChart && chart1!=mChart) {
            continue;
        }

        // TODO: have motion


        m4d::vec3 p0trans = volInvTransfMat * p0trans4D.getAsV3D(); // m4d::vec3(p0trans4D.x(1),p0trans4D.x(2),p0trans4D.x(3));
        m4d::vec3 p1trans = volInvTransfMat * p1trans4D.getAsV3D(); // m4d::vec3(p1trans4D.x(1),p1trans4D.x(2),p1trans4D.x(3));
        m4d::vec3 vtrans  = p1trans - p0trans;

        double tp0 = p0.x(0);
        double tp1 = p1.x(0);

        short entryFace,exitFace;
        if (!getTentryTexit( p0trans, p1trans,
                             tp0, tp1, time_Entry, time_Exit, entryFace, exitFace ) )  {
            continue; // kein Schnitt
        }

        tEntry = (time_Entry - tp0) / (tp1 - tp0);
        tExit  = (time_Exit  - tp0) / (tp1 - tp0);

        if (GvsRay::isIn(seg,tEntry,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,tEntry)))  {
            insecEntry.setDist( GvsRay::calcRayDist(seg,tEntry));

            m4d::vec4 dir   = p1trans4D - p0trans4D; // m4d::vec4( p1trans4D.x(0)-p0trans4D.x(0),p1trans4D.x(1)-p0trans4D.x(1), p1trans4D.x(2)-p0trans4D.x(2),p1trans4D.x(3)-p0trans4D.x(3) );
            m4d::vec4 point = p0trans4D + tEntry * dir;
            m4d::enum_coordinate_type cType = mMetric->getCoordType();
            m4d::TransCoordinates::coordTransf(m4d::enum_coordinate_cartesian,point,dir,cType,point,dir);
            insecEntry.setPoint( point );
            insecEntry.setDirection( dir );
            insecEntry.setSurface( this );
            insecEntry.setSubType( GvsSurfIntersec::Entering );

            insecEntry.setLocalPoint( p0trans + tEntry * vtrans );
            insecEntry.setLocalDirection ( vtrans );

            // insecEntry.partIndex = rayEntryFace;
            insecEntry.setRaySegNumber(seg);
            entryFound = true;
        }

        if (GvsRay::isIn(seg,tExit,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,tExit))) {
            insecExit.setDist( GvsRay::calcRayDist(seg,tExit));

            m4d::vec4 dir   = p1trans4D - p0trans4D; //m4d::vec4(p1trans4D.x(0)-p0trans4D.x(0),p1trans4D.x(1)-p0trans4D.x(1),  p1trans4D.x(2)-p0trans4D.x(2),p1trans4D.x(3)-p0trans4D.x(3));
            m4d::vec4 point = p0trans4D + tExit * dir;
            m4d::enum_coordinate_type cType = mMetric->getCoordType();
            m4d::TransCoordinates::coordTransf(m4d::enum_coordinate_cartesian,point,dir,cType,point,dir);
            insecExit.setPoint( point );
            insecExit.setDirection( dir );
            insecExit.setSurface( this );
            insecExit.setSubType( GvsSurfIntersec::Exiting );

            insecExit.setLocalPoint( p0trans + tExit * vtrans );
            insecExit.setLocalDirection( vtrans );

            //  insecExit.partIndex = rayExitFace;
            insecExit.setRaySegNumber(seg);
            exitFound = true;
        }

        if ( entryFound && exitFound ) {
            spanList.insert(insecEntry, insecExit);
            return true;
        }
    }
    return false;  // no intersection
}


void GvsSolEllipsoid::calcNormal( GvsSurfIntersec & intersec ) const {
    m4d::vec3 localNormal  = intersec.localPoint();
    m4d::vec3 globalNormal = transposeMult( volInvTransfMat, localNormal );
    //globalNormal.print();
    //intersec.setNormal( globalNormal.getNormalized() );
    intersec.setNormal( localNormal.getNormalized() );
}


void GvsSolEllipsoid ::calcTexUVParam ( GvsSurfIntersec & intersec ) const {
    static double inv2Pi = 1.0 / (2.0 * M_PI);
    static double invPi = 1.0 /  M_PI;

    double uParam, vParam;

    m4d::vec3 localNormal = intersec.localPoint();

    double phi = acos( - localNormal.x(2));
    vParam     = phi * invPi;

    if ( (fabs(localNormal.x(0) ) < GVS_EPS) && (fabs(localNormal.x(1) ) < GVS_EPS) ) {
        uParam = 0.0;
    }
    else {
        double theta = acos( localNormal.x(0)/sin(phi) ) * inv2Pi;

        if (localNormal.x(1) >= 0.0) {
            uParam = theta;
        }
        else {
            uParam = 1.0 - theta;
        }
    }
    intersec.setTexUVParam(uParam, vParam);
    //std::cerr << " uParam = " << uParam << "   vParam = " << vParam << std::endl;
}


void GvsSolEllipsoid :: calcBoundBox() {
    volBoundBox = GvsBoundBox ( volTransfMat * m4d::vec3( -1.0, -1.0, -1.0 ),
                                volTransfMat * m4d::vec3( -1.0,  1.0, -1.0 ) );

    volBoundBox.extendBoxToContain( volTransfMat * m4d::vec3( 1.0, -1.0, -1.0) );
    volBoundBox.extendBoxToContain( volTransfMat * m4d::vec3( 1.0,  1.0, -1.0) );
    volBoundBox.extendBoxToContain( volTransfMat * m4d::vec3(-1.0, -1.0,  1.0) );
    volBoundBox.extendBoxToContain( volTransfMat * m4d::vec3(-1.0,  1.0,  1.0) );
    volBoundBox.extendBoxToContain( volTransfMat * m4d::vec3( 1.0, -1.0,  1.0) );
    volBoundBox.extendBoxToContain( volTransfMat * m4d::vec3( 1.0,  1.0,  1.0) );
    //   volBoundBox.print(cerr);
}


void GvsSolEllipsoid::Print( FILE* fptr ) {
    fprintf(fptr,"SolEllipsoid {\n");
    fprintf(fptr,"\tcenter:  ");getCenter().printS(fptr);
    fprintf(fptr,"\taxlen:   ");getHalfAxisLength().printS(fptr);
    fprintf(fptr,"\tshader:  ");
    if (surfShader!=NULL) {
        surfShader->Print(fptr);
    }
    fprintf(fptr,"\tobjType: %s\n",(mObjType==local)?"local":"inCoords");
    fprintf(fptr,"\tchart:   %d\n",mChart);
    fprintf(fptr,"}\n\n");
    //TODO: print motion
}


void GvsSolEllipsoid :: calcDerivatives( GvsSurfIntersec & intersec ) const {
    //   std::cerr << "void GvsSolEllipsoid :: calcDerivatives called" << std::endl;

    double phi   = intersec.texUVParam().x(0);
    double theta = intersec.texUVParam().x(1);

    m4d::vec3 Du = volTransfMat * m4d::vec3( cos(phi) * sin(theta),
                                             sin(phi) * sin(theta),
                                             cos(theta) );

    intersec.setDerivS( Du.getNormalized() );

    m4d::vec3 Dv = intersec.normal() ^ intersec.derivS();
    intersec.setDerivT( Dv.getNormalized() );
}


bool GvsSolEllipsoid :: PtInsideEllipsoid ( const m4d::vec4& pt ) const {
    m4d::vec3 pt_3D( pt.x(1), pt.x(1), pt.x(1) );
    m4d::vec3 localPt = volInvTransfMat * pt_3D;

    return ( localPt.x(0) * localPt.x(0) +
             localPt.x(1) * localPt.x(1) +
             localPt.x(2) * localPt.x(2) <= 1.0 );
}
