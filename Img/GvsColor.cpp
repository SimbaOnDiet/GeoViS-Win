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

#include "GvsColor.h"

GvsColor RgbBlack  ( 0, 0, 0 );
GvsColor RgbWhite  ( 1, 1, 1 );
GvsColor RgbRed    ( 1, 0, 0 );
GvsColor RgbGreen  ( 0, 1, 0 );
GvsColor RgbBlue   ( 0, 0, 1 );
GvsColor RgbYellow ( 1, 1, 0 );
GvsColor RgbMagenta( 1, 0, 1 );
GvsColor RgbCyan   ( 0, 1, 1 );


GvsColor :: GvsColor ( double gray ) {
    red = green = blue = gray;
    alpha = 0.0;
    valid = true;
}

GvsColor :: GvsColor ( double r, double g, double b, double a ) {
    red   = r;
    green = g;
    blue  = b;
    alpha = a;
    valid = true;
}


GvsColor :: ~GvsColor() {
}

void GvsColor::setValid( bool valid ) {
    this->valid = valid;
}

bool GvsColor::isValid() {
    return valid;
}


double GvsColor::operator[] ( int i ) const {
    if ( i == 0 ) {
        return red;
    }
    else if ( i == 1 ) {
        return green;
    }
    else if ( i == 2 ) {
        return blue;
    }
    else if ( i == 3 ) {
        return alpha;
    }
    else  {
        std::cerr << "Error in GvsColor::operator[]: Index out of bounds." << std::endl;
        return red;
    }
}

double&
GvsColor :: operator[] ( int i )
{
    if ( i == 0 ) {
        return red;
    }
    else if ( i == 1 ) {
        return green;
    }
    else if ( i == 2 ) {
        return blue;
    }
    else if ( i == 3 ) {
        return alpha;
    }
    else
    {
        std::cerr << "Error in GvsColor::operator[]: Index out of bounds." << std::endl;
        return red;
    }
}

void
GvsColor :: operator= ( const GvsColor &col )
{
    red   = col.red;
    green = col.green;
    blue  = col.blue;
    alpha = col.alpha;

    data = col.data;
}

//-------------------------------------------------------------------------
//   friend function     m a x R g b D i s t
//-------------------------------------------------------------------------
double maxRgbDist ( const GvsColor &col0, const GvsColor &col1 )
{
    float rdist = fabs( col1.red   - col0.red   );
    float gdist = fabs( col1.green - col0.green );
    float bdist = fabs( col1.blue  - col0.blue  );

    return (rdist > gdist) ? ((rdist > bdist) ? rdist : bdist)
                           : ((gdist > bdist) ? gdist : bdist);
}

//-------------------------------------------------------------------------
//   friend function     m a n h a t t e n D i s t
//-------------------------------------------------------------------------
double manhattenDist ( const GvsColor &col0, const GvsColor &col1 )
{
    return fabs( col1.red   - col0.red   ) +
            fabs( col1.green - col0.green ) +
            fabs( col1.blue  - col0.blue  );
}

// --------------------------------------------------------------------------
//              trim
// --------------------------------------------------------------------------
GvsColor
GvsColor :: trim ()
{
    if ( red > 1.0 )
    {
        //std::cerr << "Warning: RvsColor::trim() called (red)." << std::endl;
        red = 1.0;
    }

    if ( green > 1.0 )
    {
        //std::cerr << "Warning: RvsColor::trim() called (green)." << std::endl;
        green = 1.0;
    }

    if ( blue > 1.0 )
    {
        //std::cerr << "Warning: RvsColor::trim() called (blue)." << std::endl;
        blue = 1.0;
    }

    return *this;
}


GvsColor GvsColor :: trimAlpha () {
    if ( alpha > 1.0 )   {
        //std::cerr << "Warning: RvsColor::trim() called (alpha)." << std::endl;
        alpha = 1.0;
    }
    return *this;
}


double
GvsColor :: luminance() const
{
    return 0.299 * red + 0.587 * green + 0.114 * blue;
}


// --------------------------------------------------------------------------
//              printRGB(a), all
// --------------------------------------------------------------------------
void
GvsColor :: printRGB ( std::ostream &os ) const
{
    os << red << " " << green << " " << blue << std::endl;
}

void
GvsColor :: printRGBa ( std::ostream &os ) const
{
    os << red << " " << green << " " << blue << " " << alpha << std::endl;
}

void
GvsColor :: printAll  ( std::ostream &os ) const
{
    std::cerr << "GvsColor :: printAll()...\n";
    os << red << " " << green << " " << blue << " " << alpha << "  ";
    //os << data.redShift << " ";
    //os << data.timeShift << " ";
    /*
  for (int i=0; i<5; i++)
    os << data.pos[i] << " ";

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      os << data.base[i][j] << " ";
*/
    os << std::endl;
}

void GvsColor::Print ( FILE* fptr ) const {
    fprintf(fptr,"%6.3f %6.3f %6.3f %6.3f\n",red,green,blue,alpha);
}

