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
    TransfMat.h  
 
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

/*!  
     \class  m4d::TranslateMat3D
     \brief  Translation matrix in 3 dimensions.
  
                    TranslateMat3D ( float, float, float )
                                   ( vec3 )
 
     \class  m4d::TranslateMat2D
     \brief  Translation matrix in 2 dimensions.
     
 
     \class  m4d::RotateMat3D
     \brief  Rotation matrix in 3 dimensions as 3x4 matrix
     
                    RotateMatd3D   ( vec3, float )

     \class  m4d::RotateMat3Df
     \brief  Rotation matrix in 3 dimensions as 3x3 matrix.

     \class  m4d::RotateMat2D
     \brief  Rotation matrix in 2 dimensions.
     
     \class  m4d::ScaleMat3D
     \brief  Scaling matrix in 3 dimensions.
     
                    ScaleMat3D ( float float float )
                               ( vec3 )
			       
     \class  m4d::ScaleMat2D
     \brief  Scaling matrix in 2 dimensions.			       
 
           Transformation matrices have to be Matrix<float,3,4> or 
           Matrix<float,2,3> !
*/
// -------------------------------------------------------------------------------- 

#ifndef M4D_TRANSF_MAT_H
#define M4D_TRANSF_MAT_H

#include <string>
#include <typeinfo>
#include <cassert>

#include <m4dGlobalDefs.h>

namespace m4d
{

//----------------------------------------------------------------------------
//         TranslateMat3D
//----------------------------------------------------------------------------
//
//                    (  1  0  0  tx )
//  TranslateMat3D =  (  0  1  0  ty )
//                    (  0  0  0  tz )
//
class MATH_API TranslateMat3D : public Matrix<double,3,4>
{
 public:
  TranslateMat3D ( double tx, double ty, double tz );
  TranslateMat3D ( const vec3 &translat );
};

//----------------------------------------------------------------------------
//         RotateMat3D
//----------------------------------------------------------------------------
class MATH_API RotateMat3D : public Matrix<double,3,4>
{
 public:
  RotateMat3D ( const vec3 &rotAxis, double rotAngle );
  RotateMat3D ( enum_axisID  mainAxis, double rotAngle );
};


class MATH_API RotateMat3Df : public Matrix<float,3,3>
{
 public:
  RotateMat3Df ( );
  RotateMat3Df ( const vec3f &rotAxis,  float rotAngle );
  RotateMat3Df ( enum_axisID  mainAxis, float rotAngle );
};


//----------------------------------------------------------------------------
//         ScaleMat3D
//----------------------------------------------------------------------------
//
//                ( sx  0  0  0  )
//  ScaleMat3D =  (  0 sy  0  0  )
//                (  0  0 sz  0  )
//
class MATH_API ScaleMat3D : public Matrix<double,3,4>
{
 public:
    ScaleMat3D ( double sx, double sy, double sz );
    ScaleMat3D ( const vec3 &scale);
};


//----------------------------------------------------------------------------
//         TranslateMat2D
//----------------------------------------------------------------------------
class MATH_API TranslateMat2D : public Matrix<double,2,3>
{
 public:
    TranslateMat2D ( double tx, double ty );
};

//----------------------------------------------------------------------------
//         RotateMat2D
//----------------------------------------------------------------------------
class MATH_API RotateMat2D : public Matrix<double,2,3>
{
 public:
    RotateMat2D ( double rotAngle );
    RotateMat2D ( double rotCenterX, double rotCenterY, double rotAngle );
};

//----------------------------------------------------------------------------
//         ScaleMat2D
//----------------------------------------------------------------------------
class MATH_API ScaleMat2D : public Matrix<double,2,3>
{
 public:
    ScaleMat2D ( double sx, double sy );
};

} // end namespace m4d

#endif

