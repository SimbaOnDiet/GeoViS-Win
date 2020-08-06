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
#ifndef GVS_LOCAL_TETRAD_H
#define GVS_LOCAL_TETRAD_H

#include <iostream>
#include <cassert>

#include "GvsGlobalDefs.h"
#include "Obj/GvsBase.h"
#include "Obj/GvsBoundBox.h"
#include "Obj/GvsBoundBox4D.h"


#include <metric/m4dMetric.h>  //LINK:LIB4MOTION


class GvsLocalTetrad : public GvsBase
{
public:
    GvsLocalTetrad ( );
    GvsLocalTetrad ( m4d::Metric* metric );
    GvsLocalTetrad ( const GvsLocalTetrad* lt );
    GvsLocalTetrad ( m4d::Metric* metric, const m4d::vec4 &p, const m4d::vec4 &v );
    GvsLocalTetrad ( m4d::Metric* metric,
                     const m4d::vec4 &e_0, const m4d::vec4 &e_1, const m4d::vec4 &e_2, const m4d::vec4 &e_3,
                     const m4d::vec4 &p, const bool coords=true);
    ~GvsLocalTetrad( );

    void  setLocalTetrad ( const GvsLocalTetrad &lT);

    void  setTetrad    ( const m4d::vec4 &e_0, const m4d::vec4 &e_1, const m4d::vec4 &e_2, const m4d::vec4 &e_3, bool gramSchmidt = false );
    void  setPosition  ( const m4d::vec4 &p );
    void  setPositionX ( int coord, double val );
    void  setVelocity  ( const m4d::vec4 &v );
    void  setAccel     ( const m4d::vec4 &a );
    void  setInCoords  ( const bool coords, const m4d::enum_nat_tetrad_type lft = m4d::enum_nat_tetrad_default );
    void  setLocalTime ( double tau );

    void  setLocalVel  ( const m4d::vec3 &v, const m4d::enum_nat_tetrad_type lft = m4d::enum_nat_tetrad_default );

    void  setE         ( int k, const m4d::vec4 &e_k );
    void  setTriad     ( const m4d::vec4 &e_1, const m4d::vec4 &e_2, const m4d::vec4 &e_3 );

    m4d::vec4 getE     ( int k ) const;
    void  getTetrad    ( m4d::vec4 &e_0, m4d::vec4 &e_1, m4d::vec4 &e_2, m4d::vec4 &e_3 );

    m4d::vec4 getBase  ( int k ) const;
    void  getInvTetrad ( m4d::vec4 &b_0, m4d::vec4 &b_1, m4d::vec4 &b_2, m4d::vec4 &b_3 );

    void  calcInvert   ( void );
    void  calcInvert   ( const m4d::vec4 &a0, const m4d::vec4 &a1, const m4d::vec4 &a2, const m4d::vec4 &a3,
                         m4d::vec4 &b0, m4d::vec4 &b1, m4d::vec4 &b2, m4d::vec4 &b3);

    bool  isRightHanded ( ) const;

    double getTime() const;
    double getLocalTime() const;

    m4d::vec4  getPosition() const;
    m4d::vec4  getVelocity() const;
    m4d::vec4  getAccel() const;

    bool       getInCoords ( void ) const;
    void       setLFType   ( const m4d::enum_nat_tetrad_type lftype );
    m4d::enum_nat_tetrad_type  getLFType   ( void ) const;

    void  adjustTetrad ( );  // test position and velocity and adapt the tetrad such that
                             // the base vector e0 points in the direction of motion.

    //! transform tetrad between coordinate and natural tetrad representation
    void  transformTetrad ( const bool coords, const m4d::enum_nat_tetrad_type lft = m4d::enum_nat_tetrad_default);

    // transform a point between coordinate and tetrad representation
    m4d::vec4 transToLocTetrad ( const m4d::vec4 &point ) const;
    m4d::vec4 transToCoords    ( const m4d::vec4 &point ) const;

    // three and four-velocity wrt. local tetrad
    m4d::vec4 getFourVelocity  ( const m4d::vec3 &v ) const;
    m4d::vec4 getFourVelocity  ( const m4d::vec4 &v ) const;
    m4d::vec3 getThreeVelocity ( const m4d::vec4 &u ) const;

    m4d::vec4 localToCoord     ( const m4d::vec4 &vec ) const;
    m4d::vec4 coordToLocal     ( const m4d::vec4 &vec ) const;


    void       setMetric ( m4d::Metric* metric );
    m4d::Metric* getMetric () const;

    void            setSTBoundBox  ( GvsBoundBox4D* box );
    GvsBoundBox4D*  getSTBoundBox  ( ) const;


    GvsLocalTetrad* getInterpolatedTetrad ( GvsLocalTetrad* lt0, GvsLocalTetrad* lt1, double frak );

    bool SetParam ( std::string pName, int val );
    bool SetParam ( std::string pName, double val );
    bool SetParam ( std::string pName, m4d::vec4 pt );
    bool SetParam ( std::string pName, m4d::vec3 vt );

    void printP ( ) const;
    void printS ( std::ostream &os = std::cout );

    virtual void Print  ( FILE* fptr = stderr );

    // ------ attributes ------
private:
    //! TRUE  : base std::vectors are given with respect to coordinates
    //! FALSE : base std::vectors are given with respect to a natural local tetrad defined in GvsMetric
    bool inCoords;

    //! Falls inCoords==false wird die lokale Tetrade bzgl dem in GvsMetric definierten
    //! lokalen Frame definiert. Der Typ des lokalen Frames ist durch lfType gegeben. The type of natural local tetrad, meaningful only inCoords==false
    m4d::enum_nat_tetrad_type lfType;

    //! Base vectors of the local tetrad
    m4d::vec4 e[4];        // e_(i) = e_(i)^mu \partial_mu

    //! inverse Base vectors
    m4d::vec4 base[4];

    //! local tetrad is a right-handed system (mRightHanded==true)
    bool mRightHanded;

    //! Position, velocity, and acceleration are with respect to coordinates.
    m4d::vec4 pos;
    m4d::vec4 vel;
    m4d::vec4 acc;

    //! proper time
    double mTau;

    //! local tetrad knows the metric it lives in.
    m4d::Metric* locTetradMetric;

    //! BoundingBox around the local tetrad (will be calculated in LocalCompObj)
    GvsBoundBox4D*   stBoundBox;
};

//----------------------------------------------------------------------------
//         inline getE(i), getPos
//----------------------------------------------------------------------------
inline m4d::vec4
GvsLocalTetrad :: getE( int k ) const {
    assert(k<4);
    return e[k];
}

inline m4d::vec4 GvsLocalTetrad :: getBase( int k ) const {
    assert(k<4);
    return base[k];
}

inline double GvsLocalTetrad :: getTime ( void ) const {
    return pos.x(0);
}

inline double GvsLocalTetrad :: getLocalTime ( void ) const {
    return mTau;
}

inline m4d::vec4 GvsLocalTetrad :: getPosition( void   ) const {
    return pos;
}

inline m4d::vec4 GvsLocalTetrad :: getVelocity( void   ) const {
    return vel;
}

inline m4d::vec4 GvsLocalTetrad :: getAccel( void   ) const {
    return acc;
}

inline bool GvsLocalTetrad :: isRightHanded ( void ) const {
    return mRightHanded;
}

#endif
