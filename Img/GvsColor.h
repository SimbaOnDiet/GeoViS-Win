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
#ifndef GVS_COLOR_H
#define GVS_COLOR_H

#include "GvsGlobalDefs.h"
#include <iostream>

/**
 * @brief The GvsColor class stores not only the RGB value of a pixel
 *        but also other data like redshift etc.
 */
class GvsColor
{
public:
    GvsColor ( double gray = 0.0 );
    GvsColor ( double r, double g, double b, double alpha = 1.0 );
    ~GvsColor();

    void      setValid( bool valid );
    bool      isValid();

    GvsColor  operator+ ( const GvsColor &col ) const;
    GvsColor  operator- ( const GvsColor &col ) const;
    GvsColor  operator* ( const GvsColor &col ) const;
    GvsColor  operator* ( double scale        ) const;
    GvsColor  operator/ ( double scale        ) const;
    void      operator= ( const GvsColor &col );

    void      operator+= ( const GvsColor &col );
    void      operator-= ( const GvsColor &col );
    void      operator*= ( const GvsColor &col );
    void      operator*= ( double scale        );
    void      operator/= ( double scale        );

    int       operator== ( const GvsColor &col ) const;
    int       operator!= ( const GvsColor &col ) const;
    int       operator<  ( const GvsColor &col ) const;
    int       operator<= ( const GvsColor &col ) const;
    int       operator>  ( const GvsColor &col ) const;
    int       operator>= ( const GvsColor &col ) const;

    double    operator[] ( int i ) const;
    double&   operator[] ( int i );

    GvsColor  trim       ( );
    GvsColor  trimAlpha  ( );
    double    luminance  ( ) const;

    void      printRGB  ( std::ostream &os = std::cerr ) const;
    void      printRGBa ( std::ostream &os = std::cerr ) const;
    void      printAll  ( std::ostream &os = std::cerr ) const;
    void      Print ( FILE* fptr = stderr ) const;

    // friends:
    friend  GvsColor  operator*     ( double scale, const GvsColor &col );
    friend  double    maxRgbDist    ( const GvsColor &col0, const GvsColor &col1 );
    friend  double    manhattenDist ( const GvsColor &col0, const GvsColor &col1 );


public:
    bool     valid;                     // is valid color
    double   red, green, blue, alpha;   // rgb-values
    gvsData  data;
};


/*-------------------------------------------------------------------------*/
/*                     C o n s t a n t s                                   */
/*-------------------------------------------------------------------------*/

extern GvsColor RgbBlack;
extern GvsColor RgbWhite;
extern GvsColor RgbRed;
extern GvsColor RgbGreen;
extern GvsColor RgbBlue;
extern GvsColor RgbYellow;
extern GvsColor RgbMagenta;
extern GvsColor RgbCyan;

/*-------------------------------------------------------------------------*/
/*                     I n l i n e    f u n c t i o n s                    */
/*-------------------------------------------------------------------------*/

inline GvsColor
GvsColor :: operator+ ( const GvsColor &col ) const
{
    return GvsColor ( red + col.red, green + col.green, blue + col.blue, alpha + col.alpha );
}

inline GvsColor
GvsColor :: operator- ( const GvsColor &col ) const
{
    return GvsColor ( red - col.red, green - col.green, blue - col.blue, alpha - col.alpha );
}

inline GvsColor
GvsColor :: operator* ( const GvsColor &col ) const
{
    return GvsColor( red * col.red, green * col.green, blue * col.blue, alpha * col.alpha );
}

inline GvsColor
GvsColor :: operator* ( double scale ) const
{
    return GvsColor( red * scale, green * scale, blue * scale, alpha * scale );
}

inline GvsColor
GvsColor :: operator/ ( double scale ) const
{
    return GvsColor( red / scale, green / scale, blue / scale, alpha / scale );
}

inline void
GvsColor :: operator+= ( const GvsColor &col )
{
    red += col.red;
    green += col.green;
    blue += col.blue;
    alpha += col.alpha;
}

inline void
GvsColor :: operator-= ( const GvsColor &col )
{
    red -= col.red;
    green -= col.green;
    blue -= col.blue;
    alpha -= col.alpha;
}

inline void
GvsColor :: operator*= ( const GvsColor &col )
{
    red *= col.red;
    green *= col.green;
    blue *= col.blue;
    alpha *= col.alpha;
}

inline void
GvsColor :: operator*= ( double scale )
{
    red *= scale;
    green *= scale;
    blue *= scale;
    alpha *= scale;
}

inline void
GvsColor :: operator/= ( double scale )
{
    red /= scale;
    green /= scale;
    blue /= scale;
    alpha /= scale;
}

inline int
GvsColor :: operator== ( const GvsColor &col ) const
{
    return (red == col.red) && (green == col.green) && (blue == col.blue) && (alpha == col.alpha);
}

inline int
GvsColor :: operator!= ( const GvsColor &col ) const
{
    return (red != col.red) || (green != col.green) || (blue != col.blue) || (alpha != col.alpha);
}

inline int
GvsColor :: operator < ( const GvsColor &col ) const
{
    return (red < col.red) && (green < col.green) && (blue < col.blue);
}

inline int
GvsColor :: operator <= ( const GvsColor &col ) const
{
    return (red <= col.red) && (green <= col.green) && (blue <= col.blue);
}

inline int
GvsColor :: operator > ( const GvsColor &col ) const
{
    return (red > col.red) && (green > col.green) && (blue > col.blue);
}

inline int
GvsColor :: operator >= ( const GvsColor &col ) const
{
    return (red >= col.red) && (green >= col.green) && (blue >= col.blue);
}

//friend
inline GvsColor operator* ( double scale, const GvsColor &col )
{
    return GvsColor( scale * col.red, scale * col.green, scale * col.blue, scale * col.alpha );
}


#endif
