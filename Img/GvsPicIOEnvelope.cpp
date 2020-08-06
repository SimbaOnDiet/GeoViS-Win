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

#include "GvsPicIOEnvelope.h"
#include "GvsChannelImg2D.h"

#define HAVE_LIBPPM

#ifdef HAVE_LIBTIFF
#include "GvsTiffIO.h"
#endif

#ifdef HAVE_LIBPNG
#include "GvsPngIO.h"
#endif

#ifdef HAVE_LIBPPM
#include "GvsPpmIO.h"
#endif

#include <cstring>
#include <cassert>
#include <cstdlib>


GvsPicIOEnvelope :: GvsPicIOEnvelope() {
    PicIOMgr = 0;
}


GvsPicIOEnvelope :: ~GvsPicIOEnvelope () {
    delete PicIOMgr;
}

void GvsPicIOEnvelope :: allocPicIOMgr( const char *filename )
{
    delete PicIOMgr;

#ifdef HAVE_LIBTIFF
    if ( has_extension( filename, ".tiff" ) || has_extension( filename, ".tif" ) ) {
        PicIOMgr = new GvsTiffIO;
        return;
    }
#endif // HAVE_LIBTIFF

#ifdef HAVE_LIBPNG
    if ( has_extension( filename, ".png" ) ) {
        PicIOMgr = new GvsPngIO;
        return;
    }
#endif // HAVE_LIBPNG

#ifdef HAVE_LIBPPM
    if ( has_extension( filename, ".ppm" ) ) {
        PicIOMgr = new GvsPpmIO;
        return;
    }
#endif // HAVE_LIBPPM
    fprintf(stderr,"File ending unknown: %s\n",filename);
    exit(1);
}


bool GvsPicIOEnvelope :: testImageExtension ( const std::string filename ) {
    std::string ext;
    int k = filename.rfind(".");
    ext.assign(filename,k+1,filename.size());

    bool imageExt_exists = false;
    if (ext=="ppm") {
        imageExt_exists = true;
    }
#ifdef HAVE_LIBTIFF
    else if ((ext=="tiff") || (ext=="tif")) {
        imageExt_exists = true;
    }
#endif
#ifdef HAVE_LIBPNG
    else if (ext=="png") {
        imageExt_exists = true;
    }
#endif
    return imageExt_exists;
}


bool GvsPicIOEnvelope :: readChannelImg( GvsChannelImg2D& chanImg,
                                        const char *filename ) {
    allocPicIOMgr( filename );
    chanImg.setFileName( filename );

    return PicIOMgr->readChannelImg( chanImg, filename );
}


bool GvsPicIOEnvelope :: writeChannelImg( GvsChannelImg2D& chanImg,
                                         const char *filename ) {
    allocPicIOMgr( filename );
    return PicIOMgr->writeChannelImg( chanImg, filename );
}

