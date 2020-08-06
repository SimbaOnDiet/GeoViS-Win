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

#include "Cam/Gvs2PICam.h"

#include "math/TransfMat.h"

Gvs2PICam :: Gvs2PICam() {
    viewHeading = 0.0;
    viewPitch   = 0.0;
    viewResolution = m4d::ivec2(100,100);

    AddParam("heading",gvsDT_DOUBLE);
    AddParam("pitch",gvsDT_DOUBLE);
    install();
}

Gvs2PICam :: Gvs2PICam(const double heading, const double pitch, const int res) {
    viewHeading = heading;
    viewPitch   = pitch;
    viewResolution = m4d::ivec2(res,res);

    AddParam("heading",gvsDT_DOUBLE);
    AddParam("pitch",gvsDT_DOUBLE);
    install();
}

Gvs2PICam::~Gvs2PICam() {
}


void Gvs2PICam::setHeading( const double heading ) {
    viewHeading = heading;
}

void Gvs2PICam::setPitch( const double pitch ) {
    viewPitch = pitch;
}

double Gvs2PICam::getHeading() const {
    return viewHeading;
}

double Gvs2PICam::getPitch() const {
    return viewPitch;
}


std::string Gvs2PICam :: install() {
    return std::string();
}


m4d::vec3 Gvs2PICam::GetRayDir ( const double x, const double y ) {
    double sx = 2.0*((x+0.5)/static_cast<double>(viewResolution.x(0)) - 0.5);
    double sy = 2.0*(0.5-(y+0.5)/static_cast<double>(viewResolution.x(0)));
    if (sx*sx+sy*sy>1.0) {
        return m4d::vec3();
    }

    double sz = sqrt(1.0-sx*sx-sy+sy);

    m4d::vec3 dir = m4d::vec3(sx,sy,sz);
    dir = dir.getNormalized();

    dir = m4d::RotateMat3D(m4d::axis_Z,viewHeading*DEG_TO_RAD) *
            m4d::RotateMat3D(m4d::axis_Y,-viewPitch*DEG_TO_RAD) *
            m4d::RotateMat3D(m4d::axis_Z,1.5*M_PI) * dir;
    return dir;
}


bool Gvs2PICam::SetParam ( std::string pName, double val ) {
    bool isOkay = GvsBase::SetParam(pName,val);
    if (isOkay) {
        if (pName=="heading") {
            viewHeading = val;
        } else if (pName=="pitch") {
            viewPitch = val;
        }
    }
    return isOkay;
}


void Gvs2PICam::Print( FILE* fptr ) {
    fprintf(fptr,"2PICam {\n");
    fprintf(fptr,"\theading: %f\n",viewHeading);
    fprintf(fptr,"\tpitch:   %f\n",viewPitch);
    fprintf(fptr,"\tres:     %d\n",viewResolution.x(0));
    fprintf(fptr,"\tfilter:  %s\n",GvsCamFilterNames[camFilter].c_str());
    fprintf(fptr,"}\n");
}
