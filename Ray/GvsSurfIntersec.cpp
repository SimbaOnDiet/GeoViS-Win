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

#include "Ray/GvsSurfIntersec.h"
//#include <Shader/Surface/GvsSurfaceShader.h>
#include <Obj/GvsSurface.h>
#include <Obj/GvsSceneObj.h>
#include <cmath>


GvsSurfIntersec :: GvsSurfIntersec() :
    insecSurface(NULL),
    insecLocalTetrad(NULL) {
    reset();
}

GvsSurfIntersec :: GvsSurfIntersec( const GvsSurfIntersec& s ) {
    assign(s);
}


GvsSurfIntersec :: ~GvsSurfIntersec() {
    //if ( insecSurface != NULL ){
    //  delete insecSurface;
    //  insecSurface = NULL;
    //}

    // dangerous operation, because tetrad is allocated outside this class
    // see e.g. GvsSolConvexPrim::testLocalIntersection()
    if ( insecLocalTetrad != NULL ) {
        delete insecLocalTetrad;
        insecLocalTetrad = NULL;
    }

}

GvsSurfIntersec& GvsSurfIntersec :: operator= (const GvsSurfIntersec& s ) {
    assign(s);
    return *this;
  //cerr << "Do not use GvsSurfIntersec :: operator=" << endl;
}

double GvsSurfIntersec :: dist() const {
    return insecDist;
}


GvsShader* GvsSurfIntersec :: shader () const {
    return (insecSurface == NULL) ? (GvsShader*)NULL : insecSurface->shader( *this );
}


GvsSceneObj* GvsSurfIntersec :: object () const {
    return insecSurface;
}


void GvsSurfIntersec :: reset () {
    insecDist                 = DBL_MAX;
    insecSubType              = Default;
    insecSurface              = NULL;
    insecSurfIsSelfDescribing = true;

    //...
    insecNormalIsValid        =
            insecDerivSIsValid        =
            insecDerivTIsValid        =
            insecTexUVParamAreValid   = false;

    insecLocalPointIsValid    =
            insecSurfSTParamAreValid  = false;

    insecRaySegNumber      = 0;

    partIndex              = 0;


    insecLocalInsec  = false;
    //insecLocalTetrad = NULL;
    insecLocalTime = 0.0;
}

void GvsSurfIntersec :: setSubType ( const InsecSubType& subType ) {
    insecSubType = subType;
}

void GvsSurfIntersec :: setPoint ( const m4d::vec4& point ) {
    insecPoint = point;
}

void GvsSurfIntersec :: setDirection ( const m4d::vec4 &dir ) {
    insecDirection = dir;
}

void GvsSurfIntersec :: setLocalPoint ( const m4d::vec3& lpoint ) {
    insecLocalPoint = lpoint;
    insecLocalPointIsValid = true;
}


void GvsSurfIntersec :: setNormal ( const m4d::vec3& normal ) {
    insecNormal = normal;
    insecNormNormal = normal.getNormalized();
    insecNormalIsValid = true;
}

void GvsSurfIntersec :: setDerivS ( const m4d::vec3& derivS ) {
    insecDerivS = derivS;
    insecDerivSIsValid = true;
}

void GvsSurfIntersec :: setDerivT ( const m4d::vec3& derivT ) {
    insecDerivT = derivT;
    insecDerivTIsValid = true;
}

void GvsSurfIntersec :: setTexUVParam ( const m4d::vec2& texUVParam ) {
    insecTexUVParam = texUVParam;
    insecTexUVParamAreValid = true;
}

void GvsSurfIntersec :: setTexUVParam ( double uParam, double vParam ) {
    insecTexUVParam[0] = uParam;
    insecTexUVParam[1] = vParam;
    insecTexUVParamAreValid = true;
}

void GvsSurfIntersec :: setSurfIsSelfDescribing ( bool sisd) {
    insecSurfIsSelfDescribing = sisd;
}

void GvsSurfIntersec :: setSurfSTParam ( const m4d::vec2& surfSTParam ) {
    insecSurfSTParam = surfSTParam;
    insecSurfSTParamAreValid = true;
}

void GvsSurfIntersec :: setSurfSTParam ( double sParam, double tParam ) {
    insecSurfSTParam[0] = sParam;
    insecSurfSTParam[1] = tParam;
    insecSurfSTParamAreValid = true;
}

void GvsSurfIntersec :: setDist ( double dist ) {
    insecDist = dist;
}

void GvsSurfIntersec :: setSurface ( GvsSurface* surf ) {
    insecSurface = surf;
}

void GvsSurfIntersec :: setRaySegNumber ( int segNumber ) {
    insecRaySegNumber = segNumber;
}


void GvsSurfIntersec :: invertNormal () {
    if ( ! insecNormalIsValid ) {
        insecSurface->calcNormal(*this);
    }

    insecNormal     = - insecNormal;
    insecNormNormal = - insecNormNormal;
}


GvsSurfIntersec::InsecSubType GvsSurfIntersec :: subType () const {
    return insecSubType;
}

m4d::vec4 GvsSurfIntersec :: point () const {
    return insecPoint;
}

m4d::vec4 GvsSurfIntersec::direction() const {
    return insecDirection;
}


m4d::vec3 GvsSurfIntersec :: localPoint () const {
    return insecLocalPoint;
}


m4d::vec3 GvsSurfIntersec :: normal() {
    if ( !insecNormalIsValid ) {
        insecSurface->calcNormal(*this);
    }
    return insecNormal;
}

m4d::vec3 GvsSurfIntersec :: normNormal() {
    if ( !insecNormalIsValid ) {
        insecSurface->calcNormal(*this);
    }
    return insecNormNormal;
}

m4d::vec3 GvsSurfIntersec :: derivS () {
    if ( ! insecDerivSIsValid )
        insecSurface->calcDerivatives( *this );

    return insecDerivS;
}

m4d::vec3 GvsSurfIntersec :: derivT () {
    if ( ! insecDerivTIsValid )
        insecSurface->calcDerivatives( *this );

    return insecDerivT;
}

m4d::vec2 GvsSurfIntersec :: texUVParam () {
    if ( ! insecTexUVParamAreValid )
        insecSurface->calcTexUVParam( *this );

    return insecTexUVParam;
}

m4d::vec2 GvsSurfIntersec :: surfSTParam () const {
    return insecSurfSTParam;
}

GvsSurface* GvsSurfIntersec :: surface () const {
    return insecSurface;
}

bool GvsSurfIntersec :: surfIsSelfDescribing () const {
    return insecSurfIsSelfDescribing;
}

bool GvsSurfIntersec :: surfSTParamAreValid () const {
    return insecSurfSTParamAreValid;
}

int GvsSurfIntersec :: getRaySegNumber () {
    return insecRaySegNumber;
}

void GvsSurfIntersec :: setLocalIntersec ( bool locIntersec ) {
    insecLocalInsec = locIntersec;
}

bool GvsSurfIntersec :: isLocalIntersec ( void ) const {
    return insecLocalInsec;
}


void GvsSurfIntersec :: setLocalTetrad( GvsLocalTetrad *lt ) {
    insecLocalTetrad = lt;
}

GvsLocalTetrad* GvsSurfIntersec :: getLocalTetrad() const {
    return insecLocalTetrad;
}


void GvsSurfIntersec :: setLocalDirection( const m4d::vec3& localDir ) {
    insecLocalDir = localDir;
    insecLocalDir.normalize();
}

m4d::vec3 GvsSurfIntersec :: getLocalDirection() const {
    return insecLocalDir;
}


void GvsSurfIntersec :: setLocalTime( double localTime ) {
    insecLocalTime = localTime;
}

double GvsSurfIntersec :: getLocalTime( void ) const {
    return insecLocalTime;
}


void GvsSurfIntersec :: print ( FILE* fptr ) {
    fprintf(fptr,"SurfIntersection {\n");
    fprintf(fptr,"\tdist   %f\n",insecDist);
    fprintf(fptr,"\tpoint  ");insecPoint.printS(fptr);
    fprintf(fptr,"\tdir    ");insecDirection.printS(fptr);
    fprintf(fptr,"}\n");
    /*
  os << "\tinsecDist             " << insecDist << std::endl;
  os << "\tinsecPoint            ";
  insecPoint.print(os);
  os << "\tinsecDir              ";
  insecDirection.print(os);
  os << "\tinsecSurfIsSelfDescr  " << insecSurfIsSelfDescribing << std::endl;
  os << "\tinsecRaySegNumber     " << insecRaySegNumber << std::endl;
  os << "\tinsecNormal           ";
  insecNormal.print(os);
  os << "\tinsecLocalPoint       ";
  insecLocalPoint.print(os);
  os << "\tinsecLocalDir         ";
  insecLocalDir.print(os);
  os << "}" << std::endl << std::endl;
  */
}


bool GvsSurfIntersec :: operator == ( double rayLength ) const {
    return fabs( insecDist - rayLength ) <= GVS_EPS;
}

bool GvsSurfIntersec :: operator != ( double rayLength ) const {
    return fabs( insecDist - rayLength ) > GVS_EPS;
}

bool GvsSurfIntersec :: operator < ( double rayLength ) const {
    return insecDist < (rayLength - GVS_EPS);
}

bool GvsSurfIntersec :: operator <= ( double rayLength ) const {
    return insecDist <= (rayLength + GVS_EPS);
}

bool GvsSurfIntersec :: operator > ( double rayLength ) const {
    return insecDist > (rayLength + GVS_EPS);
}

bool GvsSurfIntersec :: operator >= ( double rayLength ) const {
    return insecDist >= (rayLength - GVS_EPS);
}

bool GvsSurfIntersec :: operator == ( const GvsSurfIntersec & isec ) const {
    return operator==( isec.insecDist );
}

bool GvsSurfIntersec :: operator != ( const GvsSurfIntersec & isec ) const {
    return operator!=( isec.insecDist );
}

bool GvsSurfIntersec :: operator < ( const GvsSurfIntersec & isec ) const {
    return operator<( isec.insecDist );
}

bool GvsSurfIntersec :: operator <= ( const GvsSurfIntersec & isec ) const {
    return operator<=( isec.insecDist );
}

bool GvsSurfIntersec :: operator > ( const GvsSurfIntersec & isec ) const {
    return operator>( isec.insecDist );
}

bool GvsSurfIntersec :: operator >= ( const GvsSurfIntersec & isec ) const {
    return operator>=( isec.insecDist );
}

void GvsSurfIntersec::assign( const GvsSurfIntersec& s ) {
    this->insecDist         = s.insecDist;
    this->insecPoint        = s.insecPoint;
    this->insecDirection    = s.insecDirection;
    this->insecSubType      = s.insecSubType;
    this->insecSurface      = s.insecSurface;
    this->insecSurfIsSelfDescribing = s.insecSurfIsSelfDescribing;
    this->insecRaySegNumber = s.insecRaySegNumber;


    this->insecNormal        = s.insecNormal;
    this->insecNormNormal    = s.insecNormNormal;
    this->insecNormalIsValid = s.insecNormalIsValid;
    this->insecDerivS        = s.insecDerivS;
    this->insecDerivSIsValid = s.insecDerivSIsValid;
    this->insecDerivT        = s.insecDerivT;
    this->insecDerivTIsValid = s.insecDerivTIsValid;
    this->insecTexUVParam    = s.insecTexUVParam;
    this->insecTexUVParamAreValid = s.insecTexUVParamAreValid;

    this->insecLocalPoint        = s.insecLocalPoint;
    this->insecLocalPointIsValid = s.insecLocalPointIsValid;
    this->insecSurfSTParam       = s.insecSurfSTParam;
    this->insecSurfSTParamAreValid = s.insecSurfSTParamAreValid;
    this->insecLocalDir          = s.insecLocalDir;
    this->insecLocalTime         = s.insecLocalTime;

    this->insecLocalInsec  = s.insecLocalInsec;

    // tetrad must be new-ed, see als D'tor
    if (s.insecLocalTetrad!=NULL) {
        this->insecLocalTetrad = new GvsLocalTetrad(s.insecLocalTetrad);
    }

    this->partIndex = s.partIndex;
}
