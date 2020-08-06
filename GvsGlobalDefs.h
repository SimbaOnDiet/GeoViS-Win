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
//READING DONE
#ifndef GVS_GLOBAL_DEFS_H
#define GVS_GLOBAL_DEFS_H

#include <m4dGlobalDefs.h>

#ifdef __GNUC__
#define UNUSED_ATTR  __attribute__((unused))
#endif

//FIXED
#define UNUSED_ATTR  

#ifdef _POSIX_SOURCE
#define M_PI    3.14159265358979323846
#define M_PI_2  1.57079632679489661923
#define M_PI_4  0.78539816339744830962
#define MPITWO  6.28318530717958647692
#define M_EUL   2.71828182845904523536
#define M_EINV  0.36787944117144232160
#endif // _POSIX_SOURCE


#ifndef DBL_MAX
#define DBL_MAX 1.844674407370955616e19
#endif

#ifndef UCHAR_MAX
#define UCHAR_MAX 255
#endif


#define GVS_deg2rad  0.017453292519943295770
#define GVS_rad2deg  57.295779513082320875


#define GVS_SPEED_OF_LIGHT       2.99792458e8
#define GVS_GRAVITATIONAL_KONST  6.6742e-11
#define GVS_SOLAR_MASS           1.989e30
#define GVS_ELECTR_VACUUM_CONST  8.85418782e-12


#define GVS_ROUND(doubleval) (((doubleval) > 0.0) ? floor((doubleval)+0.5) : ceil((doubleval)-0.5))

#define GVS_TRUNC(doubleval) (((doubleval) > 0.0) ? floor((doubleval)) : ceil((doubleval)))

#define GVS_EPS    1.0e-10

#define GVS_MAX(a,b) (((a)>(b))?(a):(b))
#define GVS_MIN(a,b) (((a)<(b))?(a):(b))

#define GVS_MIN2(a,b)     GVS_MIN(a,b)
#define GVS_MAX2(a,b)     GVS_MAX(a,b)
#define GVS_MIN3(a,b,c)   ( GVS_MIN( a, GVS_MIN(b,c) ) )
#define GVS_MAX3(a,b,c)   ( GVS_MAX( a, GVS_MAX(b,c) ) )
#define GVS_MIN4(a,b,c,d) ( GVS_MIN( GVS_MIN(a,b), GVS_MIN(c,d) ) )
#define GVS_MAX4(a,b,c,d) ( GVS_MAX( GVS_MAX(a,b), GVS_MAX(c,d) ) )

#define GVS_DELTA(a,b) (((a)==(b))?(1):(0))

#define GVS_SIGN(a) ((a>0)?(1):(-1))


typedef unsigned long int   ulong;
typedef unsigned int        uint;
typedef unsigned short int  ushort;
typedef unsigned char       uchar;


enum GvsDataType {
    gvsDT_UNKNOWN = 0,
    gvsDT_INT,
    gvsDT_FLOAT,
    gvsDT_DOUBLE,
    gvsDT_VEC2,
    gvsDT_VEC3,
    gvsDT_VEC4,
    gvsDT_MAT2D,
    gvsDT_MAT3D,
    gvsDT_IVEC2,
    gvsDT_IVEC3,
    gvsDT_IVEC4,
    gvsDT_STRING,
    gvsDT_VOID
};

const std::string GvsDataTypeName[14] = {
    "gvsDT_UNKNOWN",
    "gvsDT_INT",
    "gvsDT_FLOAT",
    "gvsDT_DOUBLE",
    "gvsDT_VEC2",
    "gvsDT_VEC3",
    "gvsDT_VEC4",
    "gvsDT_MAT2D",
    "gvsDT_MAT3D",
    "gvsDT_IVEC2",
    "gvsDT_IVEC3",
    "gvsDT_IVEC4",
    "gvsDT_STRING",
    "gvsDT_VOID"
};


struct gvs_parameter
{
    GvsDataType  type;
    void*        val;
};


enum GvsCamFilter {
    gvsCamFilterRGB = 0,   // rgb image only (default)
    gvsCamFilterRGBpdz,    // rgb image + position-,direction-4-vectors + freqshift
    gvsCamFilterRGBjac     // rgb image + position-,direction-4-vectors + freqshift + Jacobi
};

const int GvsNumCamFilters = 3;
const std::string GvsCamFilterNames[GvsNumCamFilters] = {
        "FilterRGB",
        "FilterRGBpdz",
        "FilterRGBjac"
};

#define  NUM_PDZ_CHANNELS  10
#define  NUM_JAC_CHANNELS  15

typedef struct gvsData_T {
    double    objID;       // object ID of intersection object
    double    pos[4];      // position of intersection points
    double    dir[4];      // direction of light ray at intersection point (future directed)
    double    freqshift;   // gravitational frequency shift
    double    jacobi[5];   // jacobi parameters
    double    uv[2];       // uv texture coordinates
    gvsData_T() {
        objID = 0.0;
        pos[0] = pos[1] = pos[2] = pos[3] = 0.0;
        dir[0] = dir[1] = dir[2] = dir[3] = 0.0;
        freqshift = 0.0;
        jacobi[0] = jacobi[1] = jacobi[2] = jacobi[3] = jacobi[4] = 0.0;
        uv[0] = uv[1] = 0.0;
    }
} gvsData;


#endif

