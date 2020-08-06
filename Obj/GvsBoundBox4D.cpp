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

#include "Obj/GvsBoundBox4D.h"


//----------------------------------------------------------------------------
//         constructor, destructor
//----------------------------------------------------------------------------
GvsBoundBox4D::GvsBoundBox4D ( void )
{
}

GvsBoundBox4D::GvsBoundBox4D ( const m4d::vec4& lowBounds, const m4d::vec4& uppBounds,
                               bool testBounds)
{
  if ( testBounds )
  {
    for ( int i = 0; i < 4; i++ )
      if ( lowBounds[i] <= uppBounds[i] )
      {
        boxLowBounds[i] = lowBounds[i];
        boxUppBounds[i] = uppBounds[i];
      }
      else
      {
        boxLowBounds[i] = uppBounds[i];
        boxUppBounds[i] = lowBounds[i];
      }
  }
  else
  {
    for ( int i=0;i<4;i++ )
    {
      boxLowBounds[i] = lowBounds[i];
      boxUppBounds[i] = uppBounds[i];
    }
  }
}

GvsBoundBox4D :: GvsBoundBox4D ( double t,  double x,  double y,  double z,
                                 double dt, double dx, double dy, double dz,
                                 bool testBounds)
{
  if ( testBounds )
  {
    if ( dt >= 0 )
    {
      boxLowBounds.x(0) = t;
      boxUppBounds.x(0) = t + dt;
    }
    else
    {
      boxLowBounds.x(0) = t + dt;
      boxUppBounds.x(0) = t;
    }

    if ( dx >= 0 )
    {
      boxLowBounds.x(1) = x;
      boxUppBounds.x(1) = x + dx;
    }
    else
    {
      boxLowBounds.x(1) = x + dx;
      boxUppBounds.x(1) = x;
    }

    if ( dy >= 0 )
    {
      boxLowBounds.x(2) = y;
      boxUppBounds.x(2) = y + dy;
    }
    else
    {
      boxLowBounds.x(2) = y + dy;
      boxUppBounds.x(2) = y;
    }

    if ( dz >= 0 )
    {
      boxLowBounds.x(3) = z;
      boxUppBounds.x(3) = z + dz;
    }
    else
    {
      boxLowBounds.x(3) = z + dz;
      boxUppBounds.x(3) = z;
    }
  }
  else
  {
    boxLowBounds.x(0) = t;
    boxLowBounds.x(1) = x;
    boxLowBounds.x(2) = y;
    boxLowBounds.x(3) = z;
    boxUppBounds.x(0) = t + dt;
    boxUppBounds.x(1) = x + dx;
    boxUppBounds.x(2) = y + dy;
    boxUppBounds.x(3) = z + dz;
  }
}

GvsBoundBox4D::GvsBoundBox4D ( const GvsBoundBox4D& otherBox )
{
  boxLowBounds = otherBox.boxLowBounds;
  boxUppBounds = otherBox.boxUppBounds;
}

//-------------------------------------------------------------------------
//         isEmpty
//-------------------------------------------------------------------------
bool
GvsBoundBox4D :: isEmpty() const	//PROMBLEMATIC
{
  return (boxLowBounds[0] >= boxUppBounds[0]) ||
         (boxLowBounds[1] >= boxUppBounds[1]) ||
         (boxLowBounds[2] >= boxUppBounds[2]);
}


//-------------------------------------------------------------------------
//         volume
//-------------------------------------------------------------------------
double
GvsBoundBox4D :: volume() const		//PROBLEMATIC
{
  return (boxUppBounds[0] - boxLowBounds[0]) *
         (boxUppBounds[1] - boxLowBounds[1]) *
         (boxUppBounds[2] - boxLowBounds[2]);
}


//-------------------------------------------------------------------------
//         surface
//-------------------------------------------------------------------------
double
GvsBoundBox4D :: surface() const	//PROBLEMATIC
{
  m4d::vec4 size = boxUppBounds - boxLowBounds;

  return 2.0 * (size[0] * size[1] + size[1] * size[2] + size[2] * size[0]);
}

//-------------------------------------------------------------------------
//         print
//-------------------------------------------------------------------------
void GvsBoundBox4D ::  Print ( FILE* fptr ) {
    fprintf(fptr,"BoundBox4D {\n");
    fprintf(fptr,"\tlow  %7.4f %7.4f %7.4f %7.4f\n",boxLowBounds.x(0),boxLowBounds.x(1),boxLowBounds.x(2),boxLowBounds.x(3));
    fprintf(fptr,"\tupp  %7.4f %7.4f %7.4f %7.4f\n",boxUppBounds.x(0),boxUppBounds.x(1),boxUppBounds.x(2),boxUppBounds.x(3));
    fprintf(fptr,"}\n");
}
