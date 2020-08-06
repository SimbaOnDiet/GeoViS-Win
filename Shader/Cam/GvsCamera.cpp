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
//DONE READING
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

#include "GvsCamera.h"

GvsCamera::GvsCamera() : GvsBase() {
    aspectRatio = 1.0;
    camFilter = gvsCamFilterRGB;

    mRedShift = mTimeShift = mPolarisation = mAllData = mMask = false;

    // res parameter cannot be used for parallel calculation.
    // AddParam("res",gvsDT_IVEC2);
}

GvsCamera::GvsCamera( const GvsCamFilter filter ) {
    aspectRatio = 1.0;
    camFilter = filter;

    // res parameter cannot be used for parallel calculation.
    // AddParam("res",gvsDT_IVEC2);
}

GvsCamera::~GvsCamera() {
}

void GvsCamera::PixelToAngle ( const double , const double , double &, double & ) {
    fprintf(stderr,"PixelToAngle: not implemented yet\n");
}

void GvsCamera::SetResolution( const m4d::ivec2 &res ) {
    viewResolution = res;
}

m4d::ivec2 GvsCamera::GetResolution() const {
    return viewResolution;
}

void GvsCamera :: setAspectRatio(double a) {
    aspectRatio = a;
}

double GvsCamera :: getAspectRatio() const {
    return aspectRatio;
}

void GvsCamera :: setCalcEye(int b) {
    calcEye = b;
}

int GvsCamera :: getCalcEye() const {
    return calcEye;
}

void GvsCamera :: setParameter ( double param ) {
    mParam = param;
}

double GvsCamera::getParameter() {
    return mParam;
}

void GvsCamera :: setCamFilter ( GvsCamFilter filter ) {
    camFilter = filter;
    mRedShift = mTimeShift = mPolarisation = false;
}

GvsCamFilter GvsCamera :: getCamFilter ( ) const {
    return camFilter;
}

bool GvsCamera :: isRedshift() {
    return mRedShift;
}

bool GvsCamera :: isTimelapse() {
    return mTimeShift;
}

bool GvsCamera :: isPolarisation() {
    return mPolarisation;
}

bool GvsCamera :: isAllData() {
    return mAllData;
}

bool GvsCamera :: isMask() {
    return mMask;
}

bool GvsCamera::SetParam(std::string pName, m4d::ivec2 val ) {
    bool isOkay = GvsBase::SetParam(pName,val);
    return isOkay;
}

void GvsCamera::Print( FILE* fptr ) {
    fprintf(fptr,"Camera {\n");
    fprintf(fptr,"\tres  %4d x %4d\n",viewResolution.x(0),viewResolution.x(1));
    fprintf(fptr,"\tfilt %s\n",GvsCamFilterNames[camFilter].c_str());
    fprintf(fptr,"}\n");
}
