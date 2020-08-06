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
#ifndef GVS_STMOTION_H
#define GVS_STMOTION_H

#include <iostream>
#include <cstdio>
#include <cassert>
#include <deque>

#include "GvsGlobalDefs.h"
#include "Obj/GvsBase.h"
#include "Obj/GvsBoundBox4D.h"
#include "Obj/STMotion/GvsLocalTetrad.h"



enum GvsMotionType
{
  gvsMotionNoMotion, gvsMotionGeneral, gvsMotionGeodesic, gvsMotionConstVelocity, gvsMotionUnknown
};

const std::string GvsMotionTypeName[5] = {
    "NoMotion",
    "GeneralMotion",
    "GeodesicMotion",
    "ConstVelMotion",
    "UnknownMotion"
};

typedef std::deque< GvsLocalTetrad* > dequeLocT;
typedef dequeLocT::iterator      dequeLocTit;

class GvsStMotion : public GvsBase
{
  public:
    GvsStMotion();
    virtual ~GvsStMotion();

    void  setNoMotion();
    int   getNumPositions ( void   ) const;

    GvsMotionType  getMotionType ( ) const;

    m4d::vec4   getActualPosition ( double coordTime ) const;

    // access to a concrete local tetrad
    m4d::vec4   getPosition ( int k ) const;
    m4d::vec4   getVelocity ( int k ) const;
    m4d::vec4   getAccel    ( int k ) const;
    m4d::vec4   getE        ( int k, int j ) const;   // get E(j)

    m4d::vec4   getFirstPos ( void );
    m4d::vec4   getLastPos  ( void );

    void            setLocalTetrad ( GvsLocalTetrad* locT );
    GvsLocalTetrad* getLocalTetrad ( unsigned int k );

    void            setSTBoundBox  ( GvsBoundBox4D* box, int k=0 );
    GvsBoundBox4D*  getSTBoundBox  ( int k = 0 ) const;

    // delete all entries in localTetrad
    virtual void  deleteAllEntries ( );

    // get local tetrad that is closest to the requested time.
    GvsLocalTetrad*  getClosestLT ( double time, int &num );
    
    // zugehörig zur tetrade "num" (vgl. getClosestLT) gehörende
    // eigenzeit am objekt.
    // WICHTIG: keinen integrator mit schrittweitensteuerung verwenden,
    // sonst stimmt die Annahme tau=tau0+tetradNr*dTau nicht mehr!
    double           getLocalTime ( int num );

    // only used for MotionConstVelocity
    virtual bool getTransformedPolygon(const int segNr,
                                       const m4d::vec4& p0in, const m4d::vec4& p1in, m4d::vec4& p0out, m4d::vec4& p1out);


    virtual void Print    ( FILE* fptr = stderr );
    virtual void PrintAll ( FILE* fptr = stderr ) const;
    virtual void PrintToFile ( const char* filename );


    // attributes
  protected:

    int numPositions;               // number of positions
    dequeLocT     localTetrad;      // double-ended queue of all positions and orientations
    dequeLocTit   localTetradPtr;   // iterator of upper queue

    GvsMotionType mType;            // motion: nomotion, general, geodesic, unknown
    
    bool          mPastWarning;      //don't flood sterr with warnings of getClosestLT
    bool          mFutureWarning;
};


//----------------------------------------------------------------------------
//         inline functions
//----------------------------------------------------------------------------
inline int
GvsStMotion :: getNumPositions( void  ) const
{
  return numPositions;
}

inline GvsMotionType
GvsStMotion :: getMotionType ( ) const
{
  return mType;
}

#endif
