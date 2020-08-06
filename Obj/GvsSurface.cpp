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

#include "GvsSurface.h"
#include "Shader/Surface/GvsSurfaceShader.h"
#include "Ray/GvsSurfIntersec.h"


GvsSurface::GvsSurface() {
    surfShader = NULL;
}

GvsSurface::GvsSurface(GvsSurfaceShader *shader) {
    surfShader = shader;
}

GvsSurface::~GvsSurface() {
    // TODO: only set 0 or delete ??
    surfShader = NULL;
}

void GvsSurface::setShader ( GvsSurfaceShader *shader ) {
    surfShader = shader;
}

GvsSurfaceShader* GvsSurface::shader () const {
    return surfShader;
}

GvsShader* GvsSurface::shader( const GvsSurfIntersec &  ) const {
    return surfShader;
}


void GvsSurface::calcTexUVParam( GvsSurfIntersec & intersec ) const {
    if (intersec.surfSTParamAreValid()) {
        intersec.setTexUVParam( intersec.surfSTParam() );
    }
    else {
        GvsBoundBox BB = boundingBox();
        m4d::vec3 BBsize = BB.uppBounds() - BB.lowBounds();

        int i0 = BBsize.mostDominantCoord();  // i0 ist der Index fuer die
        // laengste Boxkante
        int i1 = (i0 + 1) % 3;
        int i2 = (i0 + 2) % 3;

        if (BBsize[i1] < BBsize[i2])      // i1 ist der Index fuer die zweit-
            i1 = i2;                        // laengste Boxkante

        double uParam = (intersec.point()[i0] - BB.lowBounds()[i0]) / BBsize[i0];
        double vParam = (intersec.point()[i1] - BB.lowBounds()[i1]) / BBsize[i1];

        intersec.setTexUVParam( uParam, vParam );
    }
}


//-------------------------------------------------------------------------
//   Method    GvsSurface :: c a l c D e r i v a t i v e s
//-------------------------------------------------------------------------
//   Default insbes. fuer Primitiva, die nicht explizit parametrisiert
//   sind und daher auch keine offensichtlichen Ableitungen berechnen
//   koennen.
//   Es werden die Ableitungen des Einheitskugelpunktes berechnet,
//   der intersec.normal() als Normale hat.
//-------------------------------------------------------------------------

// void GvsSurface :: calcDerivatives( GvsSurfIntersec & intersec ) const
// {
//    m4d::vec3 norm = intersec.normal();
//
//    double phi = atan2( norm.y(), norm.x() );
//    double theta = asin( norm.z() );
//
//    m4d::vec3 Du = m4d::vec3( cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta) );
//    intersec.setDerivS( Du );
//
//    m4d::vec3 Dv = intersec.normal() ^ Du;
//    intersec.setDerivT( Dv );
// }

void GvsSurface::calcDerivatives( GvsSurfIntersec & intersec ) const {
    m4d::vec3 norm = intersec.normal();
    double length = norm.normalize();

    double s = atan2( norm.x(1), norm.x(0) );
    double t = atan2( norm.x(2), sqrt( norm.x(0)*norm.x(0) + norm.x(1)*norm.x(1) ) );

    m4d::vec3 q = m4d::vec3( -cos(s)*sin(t), -sin(s)*sin(t), cos(t) );
    m4d::vec3 Ft = sqrt(length) * q;
    intersec.setDerivT( Ft );

    intersec.setDerivS( Ft ^ norm );
}

void GvsSurface::Print ( FILE *fptr ) {
    fprintf(fptr,"GvsSurface {}\n");
    // TODO
}
