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

#include "GvsPinHoleCam.h"


GvsPinHoleCam::GvsPinHoleCam() {
    viewDirection  = m4d::vec3(1.0,0.0,0.0);
    viewUpVector   = m4d::vec3(0.0,0.0,1.0);
    viewField      = m4d::vec2(60.0,48.0);
    viewResolution = m4d::ivec2(720,576);

    AddParam("dir",gvsDT_VEC3);
    AddParam("vup",gvsDT_VEC3);
    AddParam("fov",gvsDT_VEC2);
    install();
}

GvsPinHoleCam::GvsPinHoleCam(const m4d::vec3 &dir, const m4d::vec3 &vup, const m4d::vec2 &fov, const m4d::ivec2 &res)
{
    viewDirection  = dir;
    viewUpVector   = vup;
    viewField      = fov;
    viewResolution = res;

    AddParam("dir",gvsDT_VEC3);
    AddParam("vup",gvsDT_VEC3);
    AddParam("fov",gvsDT_VEC2);
    install();
}

GvsPinHoleCam::~GvsPinHoleCam() {
    DelAllParam();
}

void GvsPinHoleCam::setViewDirection ( const m4d::vec3 &dir ) {
    viewDirection = dir;
}

m4d::vec3
GvsPinHoleCam :: getViewDirection ( ) const
{
  return viewDirection;
}

void
GvsPinHoleCam :: setViewUpVector ( const m4d::vec3 &vup )
{
  viewUpVector = vup;
}

m4d::vec3
GvsPinHoleCam :: getViewUpVector( ) const
{
  return viewUpVector;
}

void
GvsPinHoleCam :: setFieldOfView (const m4d::vec2 &fov )
{
  viewField = fov;
}

m4d::vec2
GvsPinHoleCam :: getFieldOfView ( ) const
{
  return viewField;
}

std::string GvsPinHoleCam::install() {
    // std::cerr << "GvsPinHoleCam :: install()...";

    if ( viewDirection == m4d::vec3(0.0,0.0,0.0) ) {
        return std::string("Direction of view is null std::vector.");
    }

    if ( viewUpVector == m4d::vec3(0.0,0.0,0.0) ) {
        return std::string("View up std::vector is null std::vector.");
    }

    if ( viewField.x(0) <= 0.0 || viewField.x(1) <= 0.0 ) {
        return std::string("Invalid resolution.");
    }

    HorizDir = (viewDirection^viewUpVector).getNormalized();
    VertDir  = (HorizDir^viewDirection).getNormalized();

    viewHorizHalfVector = viewDirection.getNorm()
            * tan( viewField.x(0)/360.0 * M_PI )
            * HorizDir;

    viewVertHalfVector = viewDirection.getNorm()
            * tan( viewField.x(1)/360.0 * M_PI )
            * VertDir;
    //   std::cerr << "done.\n";

    // print(cerr);
    return std::string();
}


m4d::vec3 GvsPinHoleCam::GetRayDir ( const double x, const double y )  {
    double sx, sy;
    m4d::vec3 dir;

    sx = (x + 0.5) / viewResolution.x(0);
    sy = (y + 0.5) / viewResolution.x(1);
    //fprintf(stderr,"x: %f y: %f  sx: %f sy: %f  %d %d\n",x,y,sx,sy,viewResolution.x(0),viewResolution.x(1));

    dir = viewDirection
          + (2.0 * sx - 1.0)/aspectRatio * viewHorizHalfVector
          + (1.0 - 2.0 * sy) * viewVertHalfVector;
    return dir.getNormalized();
}


void GvsPinHoleCam::PixelToAngle ( const double x, const double y, double &ksi, double &chi ) {
    m4d::vec3 dir = GetRayDir(x,y);

    chi = acos((VertDir|dir));
    double cy = (HorizDir|dir);
    double cx = (viewDirection|dir);
    ksi = atan2(cy,cx);
    //  dir.printS();
}


bool GvsPinHoleCam::SetParam(std::string pName, m4d::vec2 vec ) {
    bool isOkay = GvsBase::SetParam(pName,vec);
    if (isOkay && pName=="fov") {
        viewField = vec;
        install();
    }
    return isOkay;
}

bool GvsPinHoleCam::SetParam( std::string pName, m4d::vec3 vec ) {
    bool isOkay = GvsBase::SetParam(pName,vec);
    if (isOkay) {
        if (pName=="dir") setViewDirection(vec);
        else if (pName=="vup") setViewUpVector(vec);
        install();
    }
    return isOkay;
}


void GvsPinHoleCam::Print( FILE* fptr ) {
    fprintf(fptr,"PinholeCamera {\n");
    fprintf(fptr,"\tdir    %6.3f %6.3f %6.3f\n",viewDirection.x(0),viewDirection.x(1),viewDirection.x(2));
    fprintf(fptr,"\tvup    %6.3f %6.3f\n",viewUpVector.x(0),viewUpVector.x(1));
    fprintf(fptr,"\tfov    %6.3f %6.3f\n",viewField.x(0),viewField.x(1));
    fprintf(fptr,"\tres    %4d %4d\n",viewResolution.x(0),viewResolution.x(1));
    fprintf(fptr,"\tfilter %s\n",GvsCamFilterNames[camFilter].c_str());
    fprintf(fptr,"}\n\n");
}
