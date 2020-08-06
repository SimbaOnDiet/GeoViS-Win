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

#include "Obj/STMotion/GvsStMotion.h"
#include "Obj/STMotion/GvsLocalTetrad.h"

#include <functional>
#include <algorithm>

class testTime : public std::binary_function<GvsLocalTetrad*,double,bool>
{
public:
    // less time
    bool operator()(const double val,const GvsLocalTetrad* locT ) const
    {
        return (((locT->getPosition()).x(0))>val);
    }

    // greater time
    bool operator()(const GvsLocalTetrad* locT, const double val ) const
    {
        return (((locT->getPosition()).x(0))<val);
    }
};


GvsStMotion::GvsStMotion()
{
    mPastWarning = mFutureWarning = true;
    setNoMotion();
}


GvsStMotion::~GvsStMotion()
{
    deleteAllEntries();
}


void GvsStMotion :: deleteAllEntries() {
    if (!localTetrad.empty())  {
        for (unsigned int i=0;i<localTetrad.size();i++) {
            if (localTetrad[i]!=NULL)
                delete localTetrad[i];
        }
    }

    localTetrad.clear();
    numPositions = 0;
}


void GvsStMotion :: setNoMotion() {
    if (localTetrad.size()!=0) {
        deleteAllEntries();
    }
    mType = gvsMotionNoMotion;
    numPositions = 0;
}


m4d::vec4
GvsStMotion :: getPosition ( int k  ) const
{
    assert (k < numPositions);
    return localTetrad[k]->getPosition();
}

m4d::vec4
GvsStMotion :: getVelocity ( int k  ) const
{
    assert (k < numPositions);
    return localTetrad[k]->getVelocity();
}

m4d::vec4
GvsStMotion :: getAccel    ( int k  ) const
{
    assert (k < numPositions);
    return localTetrad[k]->getAccel();
}

m4d::vec4
GvsStMotion :: getE ( int k, int j ) const
{
    assert (k < numPositions);
    return localTetrad[k]->getE(j);
}

m4d::vec4
GvsStMotion :: getFirstPos ( void )
{
    localTetradPtr = localTetrad.begin();
    return (*localTetradPtr)->getPosition();
}

m4d::vec4
GvsStMotion :: getLastPos  ( void )
{
    return (localTetrad[localTetrad.size()-1])->getPosition();
}


void
GvsStMotion :: setLocalTetrad ( GvsLocalTetrad* locT )
{
    localTetrad.push_back(locT);
    numPositions = localTetrad.size();
}

GvsLocalTetrad*
GvsStMotion :: getLocalTetrad ( unsigned int k )
{
    assert (!localTetrad.empty());
    assert (k < localTetrad.size());
    return localTetrad[k];
}


void
GvsStMotion :: setSTBoundBox  ( GvsBoundBox4D* box, int k )
{
    assert (k < numPositions);
    localTetrad[k]->setSTBoundBox(box);
}

GvsBoundBox4D*
GvsStMotion :: getSTBoundBox  ( int k ) const
{
    assert (k < numPositions);
    return localTetrad[k]->getSTBoundBox();
}


GvsLocalTetrad*
GvsStMotion :: getClosestLT ( double time, int &num ) {
    if ( (localTetrad[0]->getTime()) > time  ) {
        if ( mPastWarning ) {
            fprintf(stderr,"GvsStMotion::getClosesLT: Motion calculation backward in time is not sufficient!  %f\n",time);
            //std::cerr << "GvsStMotion::getClosesLT:  Bewegung reicht zeitlich zu wenig in die Vergangenheit! " << time << std::endl;
            mPastWarning = false;
        }
        return NULL;
    }
    else if ( (localTetrad[localTetrad.size()-1]->getTime()) < time ) {
        if ( mFutureWarning ) {
            fprintf(stderr,"GvsStMotion::getClosesLT: Motion calculation forward in time is not sufficient!  %f\n",time);
            ////std::cerr << "GvsStMotion::getClosesLT:  Bewegung reicht zeitlich zu wenig in die Zukunft! " << time << std::endl;
            mFutureWarning = false;
        }
        return NULL;
    }

    localTetradPtr = lower_bound(localTetrad.begin(),localTetrad.end(),time,testTime());
    num = localTetradPtr-localTetrad.begin();

    return (*localTetradPtr);
}


double
GvsStMotion :: getLocalTime( int num )
{
    //return mTau0 + double(num)*mDTau;
    if ( num < 0 || num > static_cast<int>(localTetrad.size())-1 )
        return 0.0;

    //  printf("%i\n",num);
    return localTetrad.at(num)->getLocalTime();
}


bool GvsStMotion :: getTransformedPolygon ( const int ,
                                            const m4d::vec4& , const m4d::vec4& , m4d::vec4& , m4d::vec4&  ) {
    std::cerr << "GvsStMotion::getTransformedPolygon() ... not implemented yet!\n";
    return false;
}


void GvsStMotion :: Print( FILE* fptr  ) {
    fprintf(fptr,"-----------------------------------------\n");
    fprintf(fptr,"              Motion:\n");
    fprintf(fptr,"-----------------------------------------\n");
    fprintf(fptr,"\tnumPositions:   %d\n",numPositions);
    fprintf(fptr,"\tMotionType:     %s\n",GvsMotionTypeName[mType].c_str());
    /*
    if (!localTetrad.empty()) {
        os << "\n\tFirst Position: " << std::endl;
        localTetrad[0]->Print();
        os << "\tLast Position: " << std::endl;
        localTetrad[localTetrad.size()-1]->Print();
    }
    */
    fprintf(fptr,"\n");
}


void GvsStMotion :: PrintAll ( FILE*  ) const {
}


void GvsStMotion :: PrintToFile ( const char* filename ) {
    m4d::vec4 pos;
    m4d::vec4 vel;
    m4d::vec4 base[4];
    double tau;

    // output is in the following format:
    //  position, velocity, e0, e1, e2, e3

    FILE* filePtr;
    filePtr = fopen(filename,"w");

    for (int i = 0; i < numPositions; i++) {
        tau = localTetrad[i]->getLocalTime();
        pos = localTetrad[i]->getPosition();
        vel = localTetrad[i]->getVelocity();
        for (int j=0; j<4; j++) {
            base[j] = localTetrad[i]->getE(j);
        }
        fprintf(filePtr,"%12.6f    %12.6f %12.6f %12.6f %12.6f    %12.6f %12.6f %12.6f %12.6f      %12.6f %12.6f %12.6f %12.6f    %12.6f %12.6f %12.6f %12.6f    %12.6f %12.6f %12.6f %12.6f    %12.6f %12.6f %12.6f %12.6f\n",
                tau,
                pos.x(0),pos.x(1),pos.x(2),pos.x(3),
                vel.x(0),vel.x(1),vel.x(2),vel.x(3),
                base[0].x(0),base[0].x(1),base[0].x(2),base[0].x(3),
                base[1].x(0),base[1].x(1),base[1].x(2),base[1].x(3),
                base[2].x(0),base[2].x(1),base[2].x(2),base[2].x(3),
                base[3].x(0),base[3].x(1),base[3].x(2),base[3].x(3) );
    }
    fclose(filePtr);
}

