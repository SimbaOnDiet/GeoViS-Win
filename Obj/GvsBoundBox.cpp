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

#include "Obj/GvsBoundBox.h"

#include "math/TransfMat.h"   //LINK:LIB4MOTION

//----------------------------------------------------------------------------
//         constructor, destructor
//----------------------------------------------------------------------------
GvsBoundBox::GvsBoundBox ( void )
{
}

GvsBoundBox::GvsBoundBox ( const m4d::vec3& lowBounds, const m4d::vec3& uppBounds,
                           bool testBounds)
{
    if ( testBounds )
    {
        for ( int i = 0; i < 3; i++ )
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
        boxLowBounds = lowBounds;
        boxUppBounds = uppBounds;
    }
}

GvsBoundBox :: GvsBoundBox ( double x,  double y,  double z,
                             double dx, double dy, double dz,
                             bool testBounds)
{
    if ( testBounds )
    {
        if ( dx >= 0 )
        {
            boxLowBounds.x(0) = x;
            boxUppBounds.x(0) = x + dx;
        }
        else
        {
            boxLowBounds.x(0) = x + dx;
            boxUppBounds.x(0) = x;
        }

        if ( dy >= 0 )
        {
            boxLowBounds.x(1) = y;
            boxUppBounds.x(1) = y + dy;
        }
        else
        {
            boxLowBounds.x(1) = y + dy;
            boxUppBounds.x(1) = y;
        }

        if ( dz >= 0 )
        {
            boxLowBounds.x(2) = z;
            boxUppBounds.x(2) = z + dz;
        }
        else
        {
            boxLowBounds.x(2) = z + dz;
            boxUppBounds.x(2) = z;
        }
    }
    else
    {
        boxLowBounds.x(0) = x;
        boxLowBounds.x(1) = y;
        boxLowBounds.x(2) = z;
        boxUppBounds.x(0) = x + dx;
        boxUppBounds.x(1) = y + dy;
        boxUppBounds.x(2) = z + dz;
    }
}

GvsBoundBox::GvsBoundBox ( const GvsBoundBox& otherBox )
{
    boxLowBounds = otherBox.boxLowBounds;
    boxUppBounds = otherBox.boxUppBounds;
}

//-------------------------------------------------------------------------
//   Method     GvsBoundBox :: i s E m p t y
//-------------------------------------------------------------------------
bool
GvsBoundBox :: isEmpty() const
{
    return (boxLowBounds[0] >= boxUppBounds[0]) ||
            (boxLowBounds[1] >= boxUppBounds[1]) ||
            (boxLowBounds[2] >= boxUppBounds[2]);
}


//-------------------------------------------------------------------------
//   Method     GvsBoundBox :: v o l u m e
//-------------------------------------------------------------------------
double
GvsBoundBox :: volume() const
{
    return (boxUppBounds[0] - boxLowBounds[0]) *
            (boxUppBounds[1] - boxLowBounds[1]) *
            (boxUppBounds[2] - boxLowBounds[2]);
}


//-------------------------------------------------------------------------
//   Method     GvsBoundBox :: s u r f a c e
//-------------------------------------------------------------------------
double
GvsBoundBox :: surface() const
{
    m4d::vec3 size = boxUppBounds - boxLowBounds;

    return 2.0 * (size[0] * size[1] + size[1] * size[2] + size[2] * size[0]);
}

//-------------------------------------------------------------------------
//     operator=
//-------------------------------------------------------------------------
bool
GvsBoundBox :: contains  ( const m4d::vec3 &point ) const
{
    return ( (point[0] >= boxLowBounds[0]) && (point[0] <= boxUppBounds[0]) &&
            (point[1] >= boxLowBounds[1]) && (point[1] <= boxUppBounds[1]) &&
            (point[2] >= boxLowBounds[2]) && (point[2] <= boxUppBounds[2]) );
}

//-------------------------------------------------------------------------
//     extendBoxToContain
//-------------------------------------------------------------------------

void
GvsBoundBox :: extendBoxToContain( const m4d::vec3& pt )
{
    for ( int i = 0; i < 3; i++ )
        if ( pt[i] < boxLowBounds[i] )
            boxLowBounds[i] = pt[i];
        else if ( pt[i] > boxUppBounds[i] )
            boxUppBounds[i] = pt[i];
}

//-------------------------------------------------------------------------
//     operator=
//-------------------------------------------------------------------------
GvsBoundBox&
GvsBoundBox :: operator= ( const GvsBoundBox &otherBox )
{
    boxLowBounds = otherBox.boxLowBounds;
    boxUppBounds = otherBox.boxUppBounds;
    return *this;
}

//-------------------------------------------------------------------------
//     operator+   ( union )
//-------------------------------------------------------------------------
GvsBoundBox
GvsBoundBox :: operator+ ( const GvsBoundBox &otherBox ) const
{
    GvsBoundBox unionBox;

    for (int i=0; i<3; i++)
    {
        unionBox.boxLowBounds[i] = (boxLowBounds[i] < otherBox.boxLowBounds[i])
                ? boxLowBounds[i]
                  : otherBox.boxLowBounds[i];
        unionBox.boxUppBounds[i] = (boxUppBounds[i] > otherBox.boxUppBounds[i])
                ? boxUppBounds[i]
                  : otherBox.boxUppBounds[i];
    }
    return unionBox;
}

//-------------------------------------------------------------------------
//     operator+   ( union )
//-------------------------------------------------------------------------
GvsBoundBox&
GvsBoundBox :: operator+= ( const GvsBoundBox &otherBox )
{
    for (int i=0; i<3; i++ )
    {
        if ( otherBox.boxLowBounds[i] < boxLowBounds[i] )
            boxLowBounds[i] = otherBox.boxLowBounds[i];
        if ( otherBox.boxUppBounds[i] > boxUppBounds[i] )
            boxUppBounds[i] = otherBox.boxUppBounds[i];
    }
    return *this;
}

//-------------------------------------------------------------------------
//     operator*   ( intersection )
//-------------------------------------------------------------------------
GvsBoundBox
GvsBoundBox :: operator* ( const GvsBoundBox &otherBox ) const
{
    GvsBoundBox intersecBox;

    for (int i=0; i<3; i++ )
    {
        intersecBox.boxLowBounds[i] = (boxLowBounds[i] > otherBox.boxLowBounds[i])
                ? boxLowBounds[i]
                  : otherBox.boxLowBounds[i];
        intersecBox.boxUppBounds[i] = (boxUppBounds[i] < otherBox.boxUppBounds[i])
                ? boxUppBounds[i]
                  : otherBox.boxUppBounds[i];
    }
    return intersecBox;
}

//-------------------------------------------------------------------------
//     operator*=   ( intersection )
//-------------------------------------------------------------------------
GvsBoundBox&
GvsBoundBox :: operator*= ( const GvsBoundBox &otherBox )
{
    for (int i=0; i<3; i++ )
    {
        if ( otherBox.boxLowBounds[i] > boxLowBounds[i] )
            boxLowBounds[i] = otherBox.boxLowBounds[i];
        if ( otherBox.boxUppBounds[i] < boxUppBounds[i] )
            boxUppBounds[i] = otherBox.boxUppBounds[i];
    }
    return *this;
}


//--------------------------------------------------------------------------------
//						SIMBA:Get entry time and exit time
//		It assumes that you have already konw the entry point and exit point
//---------------------------------------------------------------------------------
bool GvsBoundBox :: getTentryTexit ( const m4d::vec3 &p0, const m4d::vec3 &p1, double tp0, double tp1,
                                     double &time_Entry, double &time_Exit,
                                     short &entryFace, short &exitFace ) const
{
    double t1,t2;
    short coord;

    //cerr << "enter tentryexit  " << p0 << "  " << p1 << std::endl;
    double delta_tp = tp1-tp0;

    double mx = (p1.x(0) - p0.x(0)) / delta_tp;
    double my = (p1.x(1) - p0.x(1)) / delta_tp;
    double mz = (p1.x(2) - p0.x(2)) / delta_tp;

    double ax = p0.x(0) - mx*tp0;
    double ay = p0.x(1) - my*tp0;
    double az = p0.x(2) - mz*tp0;

    m4d::vec3 rayOrig(ax,ay,az);
    m4d::vec3 rayDir(mx,my,mz);

    time_Entry = -DBL_MAX;
    time_Exit  = DBL_MAX;

    for ( coord = 0; coord < 3; coord++ )
    {
        if ( fabs(rayDir[coord]) < GVS_EPS )
        {
            if ( (rayOrig[coord] < boxLowBounds[coord]) ||
                 (rayOrig[coord] > boxUppBounds[coord])    )
                return false;	//Did not enter this box at this time
        }
        else
        {
            t1 = (boxLowBounds[coord] - rayOrig[coord]) / rayDir[coord];
            t2 = (boxUppBounds[coord] - rayOrig[coord]) / rayDir[coord];

            if ( t1 < t2 )
            {
                if ( t1 > time_Entry )
                {
                    time_Entry = t1;
                    entryFace = (coord << 1);     //Remember << 1 is the same as multiply by 2, but it works faster. We do it so as to label the opposite face.
                }
                if ( t2 < time_Exit  )
                {
                    time_Exit = t2;
                    exitFace = (coord << 1) + 1;
                }
            }
            else
            {
                if ( t2 > time_Entry )
                {
                    time_Entry = t2;
                    entryFace = (coord << 1) + 1;
                }
                if ( t1 < time_Exit  )
                {
                    time_Exit = t1;
                    exitFace = (coord << 1);
                }
            }

            //if( (time_Entry > time_Exit) || (time_Exit < RVS_EPS) ) return false;
            if ( (time_Entry > time_Exit)) return false;
        }
    }
    if (tp0 < tp1)	//You know that we are doing raytracing so we do not care much about time sequences here...It works generally.
    {
        if (time_Entry > time_Exit)
        {
            double t = time_Entry;
            time_Entry = time_Exit;
            time_Exit  = t;
            short face = entryFace;
            entryFace  = exitFace;
            exitFace   = face;
        }
    }
    else
    {
        if (time_Entry < time_Exit)
        {
            double t = time_Entry;
            time_Entry = time_Exit;
            time_Exit  = t;
            short face = entryFace;
            entryFace  = exitFace;
            exitFace   = face;
        }
    }
    //cerr << "tentryexit : intersec " << std::endl;
    return true;
}


void
GvsBoundBox :: scale ( const m4d::vec3&   scaleVec )
{
    boxLowBounds *= scaleVec;
    boxUppBounds *= scaleVec;
}

void
GvsBoundBox :: translate  ( const m4d::vec3&  transVec )
{
    boxLowBounds += transVec;
    boxUppBounds += transVec;
}

void GvsBoundBox :: rotate( const m4d::vec3&   rotAxis, double rotAngle ) {
    transform( m4d::RotateMat3D(rotAxis,-rotAngle) );   //CURIOUS Why the minus sign?
}

void GvsBoundBox :: transform( const m4d::Matrix<double,3,4>& mat ) {
    //  std::cerr << "GvsBoundBox :: transform()" << std::endl;
    if ( ! mat.isIdentMat() )
    {
        m4d::vec3 low = boxLowBounds;
        m4d::vec3 upp = boxUppBounds;					//PROBLEMATIC What comment out this function
        /*
    *this = GvsBoundBox( mat * low, mat * upp );

    extendBoxToContain( mat * m4d::vec3( upp[0], low[1], low[2] ) );
    extendBoxToContain( mat * m4d::vec3( low[0], upp[1], low[2] ) );
    extendBoxToContain( mat * m4d::vec3( upp[0], upp[1], low[2] ) );
    extendBoxToContain( mat * m4d::vec3( low[0], low[1], upp[2] ) );
    extendBoxToContain( mat * m4d::vec3( upp[0], low[1], upp[2] ) );
    extendBoxToContain( mat * m4d::vec3( low[0], upp[1], upp[2] ) );
    */
    }
}


void GvsBoundBox::Print ( FILE* fptr ) {
    fprintf(fptr,"BoundBox {\n");
    fprintf(fptr,"\t low: "); boxLowBounds.printS(fptr);
    fprintf(fptr,"\t upp: "); boxUppBounds.printS(fptr);
    fprintf(fptr,"}\n");
}
