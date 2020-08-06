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
// -------------------------------------------------------------------------------- 
/*
    TransCoordinates.h  
 
  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave
 
 
   This file is part of the m4d-library.
 
   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.
 
*/

/*!  \class  m4d::TransCoordinates
     \brief  Coordinate transformations cartesian <-> spherical <-> ...
        
             each possible tansformation has been implemented in 3 ways:
               - transform point
               - transform point and single direction
               - transform point and four directions (e.g. tetrad)
             it is possible to transform single points/dirs or whole
             arrays of points/diections:
               - directly (transCoordAB)
               - to a specific coordiante system (coordTransf)
               - to cartesian coordiantes (toCartesianCoord)
  
             Code is duplicated often to minimize function calls
*/
// -------------------------------------------------------------------------------- 

#ifndef M4D_TRANS_COORDINATES_H
#define M4D_TRANS_COORDINATES_H

#include <m4dGlobalDefs.h>

namespace m4d
{

class MATH_API TransCoordinates
{
public:
  TransCoordinates();
  ~TransCoordinates();
    
  //! transform position to cartesian coordinates
  static void toCartesianCoord ( enum_coordinate_type fromCoord, 
                                 const vec4& oldPos, vec4& newPos );

  static void toCartesianCoord ( enum_coordinate_type fromCoord, 
                                 const vec4& oldPos, const vec4& oldDir,
                                 vec4& newPos, vec4& newDir );

  static void toCartesianCoord ( enum_coordinate_type fromCoord, 
                                 const vec4& oldPos,
                                 const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                 vec4& newPos,
                                 vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3 );

  //! transform position to cartesian coordinates (ray)    
  static void toCartesianCoord ( enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                                 std::vector<vec4>& newPos );    
  static void toCartesianCoord ( enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos, const std::vector<vec4>& oldDir,
                                 std::vector<vec4>& newPos, std::vector<vec4>& newDir );
  static void toCartesianCoord ( enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                                 const std::vector<vec4>& oldDir0, const std::vector<vec4>& oldDir1, 
                                 const std::vector<vec4>& oldDir2, const std::vector<vec4>& oldDir3,
                                 std::vector<vec4>& newPos,
                                 std::vector<vec4>& newDir0, std::vector<vec4>& newDir1,
                                 std::vector<vec4>& newDir2, std::vector<vec4>& newDir3 );
  
  //! coordinate transformation from 'fromCoord' to 'toCoord' coordinates
  static void coordTransf ( enum_coordinate_type fromCoord, const vec4& oldPos,
                            enum_coordinate_type toCoord,   vec4& newPos );
  static void coordTransf ( enum_coordinate_type fromCoord, const vec4& oldPos, const vec4& oldDir,
                            enum_coordinate_type toCoord,   vec4& newPos, vec4& newDir );
  static void coordTransf ( enum_coordinate_type fromCoord, const vec4& oldPos,
                            const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                            enum_coordinate_type toCoord, vec4& newPos,
                            vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3 );
  //! coordinate transformation from 'fromCoord' to 'toCoord' coordinates (ray)    
  static void coordTransf ( enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                            enum_coordinate_type toCoord,   std::vector<vec4>& newPos );  
  static void coordTransf ( enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos, const std::vector<vec4>& oldDir,
                            enum_coordinate_type toCoord, std::vector<vec4>& newPos, std::vector<vec4>& newDir );
  static void coordTransf ( enum_coordinate_type fromCoord, const std::vector<vec4>& oldPos,
                            const std::vector<vec4>& oldDir0, const std::vector<vec4>& oldDir1,
                            const std::vector<vec4>& oldDir2, const std::vector<vec4>& oldDir3,
                            enum_coordinate_type toCoord, std::vector<vec4>& newPos,
                            std::vector<vec4>& newDir0, std::vector<vec4>& newDir1,
                            std::vector<vec4>& newDir2, std::vector<vec4>& newDir3 );
  
  
  //! coordinate transformation cart->sph
  static void transCoordCartSph ( const vec4& oldPos, vec4& newPos );
  static void transCoordCartSph ( const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir );
  static void transCoordCartSph ( const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3 );
  //! coordinate transformation sph->cart
  static void transCoordSphCart ( const vec4& oldPos, vec4& newPos );
  static void transCoordSphCart ( const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir );   
  static void transCoordSphCart ( const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3 );
  //! coordinate transformation cart->cyl
  static void transCoordCartCyl ( const vec4& oldPos, vec4& newPos );
  static void transCoordCartCyl ( const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir );
  static void transCoordCartCyl ( const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3 );
  //! coordinate transformation cyl->cart
  static void transCoordCylCart ( const vec4& oldPos, vec4& newPos );
  static void transCoordCylCart ( const vec4& oldPos, const vec4& oldDir,
                                  vec4& newPos, vec4& newDir ); 
  static void transCoordCylCart ( const vec4& oldPos, const vec4& oldDir0, const vec4& oldDir1, const vec4& oldDir2, const vec4& oldDir3,
                                  vec4& newPos, vec4& newDir0, vec4& newDir1, vec4& newDir2, vec4& newDir3 );
 
  
};

} // end namespace m4d

#endif
