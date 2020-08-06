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

#include "GvsSolConvexPrim.h"
#include <Ray/GvsSurfIntersec.h>

#include "math/TransfMat.h"


GvsSolConvexPrim :: GvsSolConvexPrim(GvsSurfaceShader* shader) : GvsSolidObj(shader) {
    volTransfMat.setIdent();
    volInvTransfMat.setIdent();

    //  std::cerr << "addparam\n";
    AddParam("transform",gvsDT_MAT3D);
    mHaveSetParamTransfMat = false;

    mID = ++mObjCounter;
}


void GvsSolConvexPrim::scale ( const m4d::vec3 &scaleVec) {
    volTransfMat = m4d::ScaleMat3D ( scaleVec ) * volTransfMat;
    volInvTransfMat *= m4d::ScaleMat3D ( 1.0/scaleVec.x(0), 1.0/scaleVec.x(1), 1.0/scaleVec.x(2) );
    calcBoundBox();
}


void GvsSolConvexPrim :: translate ( const m4d::vec3 &transVec) {
    volTransfMat = m4d::TranslateMat3D ( transVec ) * volTransfMat;
    m4d::vec3  minTrans = -transVec;
    volInvTransfMat *= m4d::TranslateMat3D ( minTrans );
    calcBoundBox();
}


void GvsSolConvexPrim :: rotate ( const m4d::vec3 &rotAxis, double rotAngle) {
    volTransfMat = m4d::RotateMat3D ( rotAxis, rotAngle ) * volTransfMat;
    volInvTransfMat *= m4d::RotateMat3D ( rotAxis, -rotAngle);
    calcBoundBox();
}


void GvsSolConvexPrim :: transform ( const m4d::Matrix<double,3,4> &mat) {
    //  std::cerr << "GvsSolConvexPrim :: transform()...\n";
    if ( !mat.isIdentMat() ) {
        volTransfMat = mat * volTransfMat;
        m4d::Matrix<double,3,4> invMat = mat;
        invMat.invert();
        volInvTransfMat = volInvTransfMat * invMat;
        calcBoundBox();
    }
}


bool GvsSolConvexPrim::SetParam( std::string pName, m4d::Matrix<double,3,4> mat ) {
    bool isOkay = GvsBase::SetParam(pName,mat);
    if (getLowCase(pName)=="transform")   {
        setParamTransfMat = mat;
        setParamInvTransfMat = mat;
        setParamInvTransfMat.invert();
        mHaveSetParamTransfMat = true;
    }
    return isOkay;
}


bool GvsSolConvexPrim :: haveSetParamTransfMat () const {
    return mHaveSetParamTransfMat;
}


GvsBoundBox GvsSolConvexPrim :: boundingBox ( ) const {
    return volBoundBox;
}

/*
bool
GvsSolConvexPrim :: getTentryTexit  ( const m4d::vec3& p0, const m4d::vec3& p1,
                                      double tp0, double tp1,
                                      double& time_Entry, double& time_Exit,
                                      short &entryFace, short &exitFace) const
{
  std::cerr << "GvsSolConvexPrim :: getTentryTexit... not implemented yet.\n";
  return false;
}
*/

bool
GvsSolConvexPrim :: getTentryTexit  ( const m4d::vec3& , const m4d::vec3& ,
                                      double , double ,
                                      double& , double& ,
                                      short &, short &) const
{
    std::cerr << "GvsSolConvexPrim :: getTentryTexit... not implemented yet.\n";
    return false;
}


/**
 * @param ray
 * @return
 */
bool GvsSolConvexPrim::testIntersection ( GvsRay &ray ) {
    assert( mMetric != NULL );

    double time_Entry, time_Exit;
    double tEntry, tExit;
    short entryFace, exitFace;
    int chart0,chart1;

    m4d::enum_coordinate_type  coords = mMetric->getCoordType();

    int maxSeg   = ray.getNumPoints()-2;
    int startSeg = int(0);
    int endSeg   = maxSeg;

    // --- loop over all segments of the ray
    for (int seg = startSeg; seg < endSeg; seg++) {
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

        m4d::vec3 p0trans = volInvTransfMat * p0trans4D.getAsV3D();
        m4d::vec3 p1trans = volInvTransfMat * p1trans4D.getAsV3D();

        if (mHaveSetParamTransfMat)  {
            // ???
            std::cerr << "uuups\n";
            p0trans = setParamInvTransfMat * p0trans;
            p1trans = setParamInvTransfMat * p1trans;
        }
        m4d::vec3 vtrans = p1trans - p0trans;

        double tp0 = p0trans4D.x(0);
        double tp1 = p1trans4D.x(0);
        if (!getTentryTexit( p0trans,p1trans, tp0,tp1, time_Entry, time_Exit, entryFace, exitFace) ) {
            continue;
        }

        tEntry = (time_Entry - tp0) / (tp1 - tp0);
        tExit  = (time_Exit  - tp0) / (tp1 - tp0);

        if (GvsRay::isIn(seg,tEntry,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,tEntry))) {
            //std::cerr << "Schnitt tEntry: " << tEntry << std::endl;

            GvsSurfIntersec surfIntersec;
            surfIntersec.setDist ( GvsRay::calcRayDist(seg,tEntry) );
            surfIntersec.setSurface ( this );

            // global intersection point and direction in proper metric coordinates
            m4d::vec4 point = p0trans4D + tEntry * vtrans4D;
            m4d::vec4 dir   = vtrans4D;
            m4d::enum_coordinate_type cType = mMetric->getCoordType();
            m4d::TransCoordinates::coordTransf(m4d::enum_coordinate_cartesian,point,dir,cType,point,dir);
            surfIntersec.setPoint( point );
            surfIntersec.setDirection( dir );

            // local intersection point in standard object system
            surfIntersec.setLocalPoint( p0trans + tEntry * vtrans );
            surfIntersec.setLocalDirection( vtrans );

            surfIntersec.partIndex = entryFace;
            surfIntersec.setRaySegNumber(seg);
            return ray.store( surfIntersec );
        }
        else if (GvsRay::isIn(seg,tExit,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,tExit))) {
            //std::cerr << "Exit-Schnitt: " << tExit << std::endl;

            GvsSurfIntersec surfIntersec;
            surfIntersec.setDist ( GvsRay::calcRayDist(seg,tExit) );
            surfIntersec.setSurface ( this );

            // global intersection point and direction in proper metric coordinates
            m4d::vec4 point = p0trans4D + tEntry * vtrans4D;
            m4d::vec4 dir   = vtrans4D;
            m4d::enum_coordinate_type cType = mMetric->getCoordType();
            m4d::TransCoordinates::coordTransf(m4d::enum_coordinate_cartesian,point,dir,cType,point,dir);
            surfIntersec.setPoint( point );
            surfIntersec.setDirection( dir );

            // local intersection point in standard object system
            surfIntersec.setLocalPoint( p0trans + tExit * vtrans );
            surfIntersec.setLocalDirection( vtrans );

            surfIntersec.partIndex = exitFace;
            surfIntersec.setRaySegNumber(seg);
            return ray.store( surfIntersec );

        }
    }
    return false;
}


bool GvsSolConvexPrim::testLocalIntersection( GvsRay &ray, const int seg,
                                              GvsLocalTetrad* lt0, GvsLocalTetrad* lt1, const m4d::vec4 p0, const m4d::vec4 p1,
                                              double& tau0, double& tau1)
{
    // std::cerr << "GvsSolConvexPrim :: testLocalIntersection\n";

    assert (mObjType==local);

    double time_Entry, time_Exit;
    double tEntry, tExit;
    short entryFace, exitFace;


    double tp0 = p0.x(0);
    double tp1 = p1.x(0);

    m4d::vec3 p0trans,p1trans;

    int maxSeg   = ray.getNumPoints()-2;

    if (stMotion != NULL)
    {
        // std::cerr << "GvsSolConvexPrim :: have motion...\n";
        if (stMotion->getMotionType() == gvsMotionConstVelocity )
        {
            m4d::vec4 p0trans4D;
            m4d::vec4 p1trans4D;
            stMotion->getTransformedPolygon(seg,p0,p1, p0trans4D, p1trans4D);
            p0trans = volInvTransfMat * m4d::vec3(p0trans4D.x(1),p0trans4D.x(2),p0trans4D.x(3));
            p1trans = volInvTransfMat * m4d::vec3(p1trans4D.x(1),p1trans4D.x(2),p1trans4D.x(3));
        }
        else {
            std::cerr << "GvsSolConvexPrim :: testLocalIntersection()...Motion not allowed...\n";
        }
    }
    else {
        p0trans = volInvTransfMat * m4d::vec3(p0.x(1),p0.x(2),p0.x(3));
        p1trans = volInvTransfMat * m4d::vec3(p1.x(1),p1.x(2),p1.x(3));
    }


    if (mHaveSetParamTransfMat) {
        p0trans = setParamInvTransfMat * p0trans;
        p1trans = setParamInvTransfMat * p1trans;
    }

    m4d::vec3 vtrans = p1trans - p0trans;
    if (!getTentryTexit(p0trans,p1trans,
                        tp0,tp1, time_Entry, time_Exit,
                        entryFace, exitFace))
    {
        return false;
    }

    tEntry = (time_Entry-tp0)/(tp1-tp0);
    tExit  = (time_Exit-tp0)/(tp1-tp0);

    if (GvsRay::isIn(seg,tEntry,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,tEntry)))
    {
        // std::cerr << "Schnitt tEntry: " << tEntry << std::endl;

        GvsSurfIntersec surfIntersec;
        surfIntersec.setDist( GvsRay::calcRayDist(seg,tEntry) );

        m4d::vec4 vp0 = ray.getPoint(seg);
        m4d::vec4 vp1 = ray.getPoint(seg+1);

        m4d::vec4 vtrans4D = m4d::vec4(vp1.x(0)-vp0.x(0),vp1.x(1)-vp0.x(1),
                                       vp1.x(2)-vp0.x(2),vp1.x(3)-vp0.x(3));
        m4d::vec4 point = m4d::vec4(vp0.x(0),vp0.x(1),vp0.x(2),vp0.x(3)) + tEntry * vtrans4D;

        //  stMetric->coordTransf( m4d::enum_coordinate_cartesian,point,stMetric->getCoordType(),point);
        surfIntersec.setPoint   ( point );
        surfIntersec.setSurface ( this );
        surfIntersec.setDirection ( vtrans4D );

        surfIntersec.setLocalIntersec ( true );
        double i;
        double frak = modf(surfIntersec.dist(),&i) ;
        GvsLocalTetrad* lt = lt0->getInterpolatedTetrad(lt0,lt1,frak);
        surfIntersec.setLocalTetrad( lt );
        surfIntersec.setLocalTime( (1.0-frak)*tau0 + frak*tau1 );
        surfIntersec.setLocalPoint( p0trans + tEntry*vtrans );
        surfIntersec.setLocalDirection ( vtrans );
        surfIntersec.partIndex = entryFace;

        surfIntersec.setRaySegNumber(seg);
        return ray.store( surfIntersec );
    }
    else if (GvsRay::isIn(seg,tExit,maxSeg) && ray.isValidSurfIntersec( GvsRay::calcRayDist(seg,tExit)))
    {
        //std::cerr << "Schnitt tExit: " << tExit << std::endl;

        GvsSurfIntersec surfIntersec;
        surfIntersec.setDist( GvsRay::calcRayDist(seg,tExit) );

        m4d::vec4 vp0 = ray.getPoint(seg);
        m4d::vec4 vp1 = ray.getPoint(seg+1);

        m4d::vec4 vtrans4D = m4d::vec4(vp1.x(0)-vp0.x(0),vp1.x(1)-vp0.x(1),
                                       vp1.x(2)-vp0.x(2),vp1.x(3)-vp0.x(3));
        m4d::vec4 point = m4d::vec4(vp0.x(0),vp0.x(1),vp0.x(2),vp0.x(3)) + tExit * vtrans4D;

        //stMetric->coordTransf( m4d::enum_coordinate_cartesian,point,stMetric->getCoordType(),point);
        surfIntersec.setPoint   ( point );
        surfIntersec.setSurface ( this );
        surfIntersec.setDirection ( vtrans4D );

        surfIntersec.setLocalIntersec ( true );
        double i;
        double frak = modf(surfIntersec.dist(),&i) ;
        GvsLocalTetrad* lt = lt0->getInterpolatedTetrad(lt0,lt1,frak);
        surfIntersec.setLocalTetrad(lt);
        surfIntersec.setLocalTime( (1.0-frak)*tau0 + frak*tau1 );
        surfIntersec.setLocalPoint( p0trans + tExit*vtrans );
        surfIntersec.setLocalDirection ( vtrans );
        surfIntersec.partIndex = exitFace;

        surfIntersec.setRaySegNumber(seg);
        return ray.store( surfIntersec );
    }
    return false;
}

