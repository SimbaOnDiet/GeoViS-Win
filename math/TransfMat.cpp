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
// -------------------------------------------------------------------------------
/*
    TransfMat.cpp

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
// -------------------------------------------------------------------------------

#include "TransfMat.h"

#ifdef _WIN32
  #ifndef __GNUC__
    #pragma warning (disable: 4244 )
  #endif
#endif


namespace m4d
{

//----------------------------------------------------------------------------
//     TranslateMat3D:    constructor, destructor
//----------------------------------------------------------------------------
TranslateMat3D :: TranslateMat3D ( double tx, double ty, double tz )
{
    setIdent();

    mat[0][3] = tx;
    mat[1][3] = ty;
    mat[2][3] = tz;

    matIsUnit = (tx == 0.0) && (ty == 0.0) && (tz == 0.0);
}

TranslateMat3D :: TranslateMat3D ( const vec3& translat )
{
   setIdent();

   mat[0][3] = translat[0];
   mat[1][3] = translat[1];
   mat[2][3] = translat[2];

   matIsUnit = (mat[0][3] == 0.0) && (mat[1][3] == 0.0) && (mat[2][3] == 0.0);
}

//----------------------------------------------------------------------------
//     RotateMat3D:    constructor, destructor
//----------------------------------------------------------------------------
//   algorithm by:  Michael E. Pique, Rotation Tools, pp. 465 - 469
//             in:  A. Glassner (ed.), Graphics Gems, Academic Press 1990
//-------------------------------------------------------------------------
RotateMat3D :: RotateMat3D ( const vec3 &rotAxis, double rotAngle )
{
    vec3 axis = rotAxis.getNormalized();

    double sinAngle = sin( rotAngle );
    double cosAngle = cos( rotAngle );
    double one_cosAngle = 1.0 - cosAngle;

    setIdent();

    mat[0][0] = one_cosAngle * axis.x(0) * axis.x(0) + cosAngle;
    mat[0][1] = one_cosAngle * axis.x(0) * axis.x(1) - sinAngle * axis.x(2);
    mat[0][2] = one_cosAngle * axis.x(0) * axis.x(2) + sinAngle * axis.x(1);

    mat[1][0] = one_cosAngle * axis.x(1) * axis.x(0) + sinAngle * axis.x(2);
    mat[1][1] = one_cosAngle * axis.x(1) * axis.x(1) + cosAngle;
    mat[1][2] = one_cosAngle * axis.x(1) * axis.x(2) - sinAngle * axis.x(0);

    mat[2][0] = one_cosAngle * axis.x(2) * axis.x(0) - sinAngle * axis.x(1);
    mat[2][1] = one_cosAngle * axis.x(2) * axis.x(1) + sinAngle * axis.x(0);
    mat[2][2] = one_cosAngle * axis.x(2) * axis.x(2) + cosAngle;

    matIsUnit = (mat[0][0]==1.0) && (mat[0][1]==0.0) && (mat[0][2]==0.0) &&
                (mat[1][0]==0.0) && (mat[1][1]==1.0) && (mat[1][2]==0.0) &&
                (mat[2][0]==0.0) && (mat[2][1]==0.0) && (mat[2][2]==1.0);
}

RotateMat3D :: RotateMat3D ( enum_axisID mainAxis, double rotAngle )
{
   setIdent();

   switch( mainAxis )
   {
      case axis_X:  mat[1][1] =  mat[2][2] = cos( rotAngle ) ;
                    mat[1][2] = - ( mat[2][1] = sin( rotAngle ) ) ;
                    matIsUnit = (mat[1][1] == 1.0) && (mat[1][2] == 0.0);
                    break;
      case axis_Y:  mat[0][0] = mat[2][2] = cos( rotAngle ) ;
                    mat[2][0] = - ( mat[0][2] = sin( rotAngle ) ) ;
                    matIsUnit = (mat[0][0] == 1.0) && (mat[2][0] == 0.0);
                    break;
      case axis_Z:  mat[0][0] = mat[1][1] = cos( rotAngle ) ;
                    mat[0][1] = - ( mat[1][0] = sin( rotAngle ) ) ;
                    matIsUnit = (mat[0][0] == 1.0) && (mat[0][1] == 0.0);
                    break;
   }
}



//----------------------------------------------------------------------------
//     RotateMat3Df:    constructor, destructor
//----------------------------------------------------------------------------
//   algorithm by:  Michael E. Pique, Rotation Tools, pp. 465 - 469
//             in:  A. Glassner (ed.), Graphics Gems, Academic Press 1990
//-------------------------------------------------------------------------
RotateMat3Df :: RotateMat3Df ( )
{
  setIdent();
}

RotateMat3Df :: RotateMat3Df ( const vec3f &rotAxis, float rotAngle )
{
    vec3f axis = rotAxis.getNormalized();

    float sinAngle = sinf( rotAngle );
    float cosAngle = cosf( rotAngle );
    float one_cosAngle = 1.0 - cosAngle;

    setIdent();

    mat[0][0] = one_cosAngle * axis.x(0) * axis.x(0) + cosAngle;
    mat[0][1] = one_cosAngle * axis.x(0) * axis.x(1) - sinAngle * axis.x(2);
    mat[0][2] = one_cosAngle * axis.x(0) * axis.x(2) + sinAngle * axis.x(1);

    mat[1][0] = one_cosAngle * axis.x(1) * axis.x(0) + sinAngle * axis.x(2);
    mat[1][1] = one_cosAngle * axis.x(1) * axis.x(1) + cosAngle;
    mat[1][2] = one_cosAngle * axis.x(1) * axis.x(2) - sinAngle * axis.x(0);

    mat[2][0] = one_cosAngle * axis.x(2) * axis.x(0) - sinAngle * axis.x(1);
    mat[2][1] = one_cosAngle * axis.x(2) * axis.x(1) + sinAngle * axis.x(0);
    mat[2][2] = one_cosAngle * axis.x(2) * axis.x(2) + cosAngle;

    matIsUnit = (mat[0][0]==1.0) && (mat[0][1]==0.0) && (mat[0][2]==0.0) &&
    (mat[1][0]==0.0) && (mat[1][1]==1.0) && (mat[1][2]==0.0) &&
    (mat[2][0]==0.0) && (mat[2][1]==0.0) && (mat[2][2]==1.0);
}

RotateMat3Df :: RotateMat3Df ( enum_axisID mainAxis, float rotAngle )
{
   setIdent();

   switch( mainAxis )
   {
      case axis_X:  mat[1][1] =  mat[2][2] = cos( rotAngle ) ;
                    mat[1][2] = - ( mat[2][1] = sin( rotAngle ) ) ;
                    matIsUnit = (mat[1][1] == 1.0) && (mat[1][2] == 0.0);
                    break;
      case axis_Y:  mat[0][0] = mat[2][2] = cos( rotAngle ) ;
                    mat[2][0] = - ( mat[0][2] = sin( rotAngle ) ) ;
                    matIsUnit = (mat[0][0] == 1.0) && (mat[2][0] == 0.0);
                    break;
      case axis_Z:  mat[0][0] = mat[1][1] = cos( rotAngle ) ;
                    mat[0][1] = - ( mat[1][0] = sin( rotAngle ) ) ;
                    matIsUnit = (mat[0][0] == 1.0) && (mat[0][1] == 0.0);
                    break;
   }
}

//----------------------------------------------------------------------------
//       ScaleMat3D:   constructor, destructor
//----------------------------------------------------------------------------
ScaleMat3D :: ScaleMat3D ( double sx, double sy, double sz )
{
   setIdent();

   mat[0][0] = sx;
   mat[1][1] = sy;
   mat[2][2] = sz;

   matIsUnit = (sx == 1.0) && (sy == 1.0) && (sz == 1.0);
}

ScaleMat3D :: ScaleMat3D ( const vec3 &scale )
{
   setIdent();

   mat[0][0] = scale[0];
   mat[1][1] = scale[1];
   mat[2][2] = scale[2];

   matIsUnit = (mat[0][0] == 1.0) && (mat[1][1] == 1.0) && (mat[2][2] == 1.0);
}



//----------------------------------------------------------------------------
//     TranslateMat2D:    constructor, destructor
//----------------------------------------------------------------------------
TranslateMat2D :: TranslateMat2D ( double tx, double ty )
{
    setIdent();

    mat[0][2] = tx;
    mat[1][2] = ty;

    matIsUnit = (tx == 0.0) && (ty == 0.0);
}

//-------------------------------------------------------------------------
//    RotateMat2D :: constructor, destructor
//-------------------------------------------------------------------------

RotateMat2D :: RotateMat2D ( double rotAngle )
{
   setIdent();

   mat[0][0] = mat[1][1] = cos( rotAngle ) ;
   mat[0][1] = - ( mat[1][0] = sin( rotAngle ) ) ;     // ?????
}

RotateMat2D :: RotateMat2D( double rotCenterX, double rotCenterY, double rotAngle )
{
    Matrix<double,2,3>  trt = TranslateMat2D ( rotCenterX, rotCenterY ) *
                              RotateMat2D ( rotAngle ) *
                              TranslateMat2D ( -rotCenterX, -rotCenterY );

    int  i ;
    int  j ;
    for(  i = 0; i < 2; i++ )
        for(  j = 0; j < 3; j++ )
            mat[i][j] = trt.getCoeff(i,j);

}


//-------------------------------------------------------------------------
//   ScaleMat2D :: constructor, destructor
//-------------------------------------------------------------------------

ScaleMat2D :: ScaleMat2D( double sx, double sy )
{
   setIdent();

   mat[0][0] = sx;
   mat[1][1] = sy;
}

} // end namespace m4d

