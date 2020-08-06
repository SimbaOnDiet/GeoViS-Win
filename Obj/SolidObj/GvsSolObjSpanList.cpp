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

#include "GvsSolObjSpanList.h"
#include <stdlib.h>    // wg. NULL !
#include <float.h>
#include <assert.h>


GvsSolObjSpanList::GvsSolObjSpanList()
{
    spanLo = new GvsSurfIntersec[ DefaultListLength ];
    assert( spanLo != NULL );

    spanHi = new GvsSurfIntersec[ DefaultListLength ];
    assert( spanHi != NULL );

    spanListLength = DefaultListLength;
    spanPtr = 0;
}


GvsSolObjSpanList::GvsSolObjSpanList( int listLength )
{
    spanListLength = (listLength > 0) ? listLength : DefaultListLength;

    spanLo = new GvsSurfIntersec[ spanListLength ];
    assert( spanLo != NULL );

    spanHi = new GvsSurfIntersec[ spanListLength ];
    assert( spanHi != NULL );

    spanPtr = 0;
}


GvsSolObjSpanList::GvsSolObjSpanList( const GvsSolObjSpanList& sl  )
{
    spanListLength = sl.spanListLength;

    spanLo = new GvsSurfIntersec[ spanListLength ];
    assert( spanLo != NULL );

    spanHi = new GvsSurfIntersec[ spanListLength ];
    assert( spanHi != NULL );

    for ( spanPtr = 0; spanPtr < sl.spanPtr; spanPtr++ ) {
        spanLo[spanPtr] = sl.spanLo[spanPtr];
        spanHi[spanPtr] = sl.spanHi[spanPtr];
    }
}


GvsSolObjSpanList::~GvsSolObjSpanList() {
    delete [] spanLo;
    delete [] spanHi;
    spanLo = spanHi = NULL;
    spanPtr = spanListLength = 0;
}


bool GvsSolObjSpanList::insert( GvsSurfIntersec &s, GvsSurfIntersec &e )
{
    GvsSurfIntersec min, max;

    if ( s < e ) {
        min = s;
        max = e;
    }
    else {
        min = e;
        max = s;
    }

    if ( spanPtr == spanListLength ) {
        if ( !increaseNumSpans(1) ) {
            return false;
        }
    }

    for( int i = 0; i < spanPtr; i++ ) {
        if ( min < spanLo[i] ) {
            for ( int j = spanPtr; j > i; j-- ) {
                spanLo[j] = spanLo[j-1];
                spanHi[j] = spanHi[j-1];
            }
            spanPtr++;
            spanLo[i] = min;
            spanHi[i] = max;
            return true;
        }
    }

    spanLo[spanPtr] = min;
    spanHi[spanPtr] = max;
    spanPtr++;
    return true;
}


int GvsSolObjSpanList::deleteFirstSpan() {
    if ( spanPtr == 0 ) {
        return -1;
    }

    for ( int i = 0; i < spanPtr-1; i++ ) {
        spanLo[i] = spanLo[i+1];
        spanHi[i] = spanHi[i+1];
    }

    return --spanPtr;
}


int GvsSolObjSpanList::deleteSpan( int i ) {
    if ( i < 0 || i >= spanPtr ) {
        return -1;
    }

    for ( int j = i; j < spanPtr-1; j++ ) {
        spanLo[j] = spanLo[j+1];
        spanHi[j] = spanHi[j+1];
    }
    return --spanPtr;
}


bool GvsSolObjSpanList::increaseNumSpans( int addnum )
{
    GvsSurfIntersec *newSpanLo = new GvsSurfIntersec[ spanListLength + addnum ];
    if ( newSpanLo == NULL ) return false;


    GvsSurfIntersec *newSpanHi = new GvsSurfIntersec[ spanListLength + addnum ];
    if ( newSpanHi == NULL ) return false;

    spanListLength += addnum;

    for ( int i = 0; i < spanPtr; i++ ) {
        newSpanLo[i] = spanLo[i];
        newSpanHi[i] = spanHi[i];
    }

    delete [] spanLo;    spanLo = newSpanLo;
    delete [] spanHi;    spanHi = newSpanHi;
    return true;
}


bool GvsSolObjSpanList::getFirstValidSpan(GvsSurfIntersec &spanStart,
                                          GvsSurfIntersec &spanEnd    ) const
{
    for ( int i = 0; i < spanPtr; i++ ) {
        if ( spanLo[i] > GVS_EPS ) {
            spanStart = spanLo[i];
            spanEnd   = spanHi[i];
            return true;
        }
    }
    return false;
}

GvsSolObjSpanList& GvsSolObjSpanList::operator=( const GvsSolObjSpanList& sl ) {
    delete [] spanLo;
    delete [] spanHi;
    spanListLength = sl.spanListLength;

    spanLo = new GvsSurfIntersec[ spanListLength ];
    assert( spanLo != NULL );

    spanHi = new GvsSurfIntersec[ spanListLength ];
    assert( spanHi != NULL );

    for ( spanPtr = 0; spanPtr < sl.spanPtr; spanPtr++ ) {
        spanLo[spanPtr] = sl.spanLo[spanPtr];
        spanHi[spanPtr] = sl.spanHi[spanPtr];
    }
    return *this;
}


GvsSolObjSpanList GvsSolObjSpanList::operator+( const GvsSolObjSpanList& rsl ) const {
    GvsSolObjSpanList slUnion( spanPtr + rsl.spanPtr );

    register int i = 0, j = 0, k = 0;
    while( i < spanPtr || j < rsl.spanPtr ) {
        if ( i < spanPtr && (j >= rsl.spanPtr || spanLo[i] < rsl.spanLo[j]) ) {
            slUnion.spanLo[k] = spanLo[i];
            slUnion.spanHi[k] = spanHi[i];
            i++;
        }
        else {
            slUnion.spanLo[k] = rsl.spanLo[j];
            slUnion.spanHi[k] = rsl.spanHi[j];
            j++;
        }

        bool spanEndFound = false;
        while ( !spanEndFound ) {
            while( i < spanPtr && spanLo[i] <= slUnion.spanHi[k] ) {
                if ( spanHi[i] > slUnion.spanHi[k] ) {
                    slUnion.spanHi[k] = spanHi[i];
                }
                i++;
            }

            spanEndFound = true;
            while( j < rsl.spanPtr && rsl.spanLo[j] <= slUnion.spanHi[k] ) {
                if ( rsl.spanHi[j] > slUnion.spanHi[k] ) {
                    slUnion.spanHi[k] = rsl.spanHi[j];
                    spanEndFound = false;
                }
                j++;
            }
        }
        k++;
    }
    slUnion.spanPtr = k;
    return slUnion;
}


GvsSolObjSpanList GvsSolObjSpanList::operator*( const GvsSolObjSpanList& rsl ) const {
    GvsSolObjSpanList isect( spanPtr + rsl.spanPtr );

    register int i = 0, j = 0, k = 0;
    while( i < spanPtr && j < rsl.spanPtr ) {
        if ( spanLo[i] <= rsl.spanLo[j] ) {
            if ( spanHi[i] < rsl.spanLo[j] ) {
                i++;
            }
            else if ( spanHi[i] < rsl.spanHi[j] ) {
                isect.spanLo[k] = rsl.spanLo[j];
                isect.spanHi[k] = spanHi[i];
                i++;
                k++;
            }
            else {
                isect.spanLo[k] = rsl.spanLo[j];
                isect.spanHi[k] = rsl.spanHi[j];
                j++;
                k++;
            }
        }
        else {
            if ( spanLo[i] > rsl.spanHi[j] ) {
                j++;
            }
            else if ( spanHi[i] > rsl.spanHi[j] ) {
                isect.spanLo[k] = spanLo[i];
                isect.spanHi[k] = rsl.spanHi[j];
                j++;
                k++;
            }
            else {
                isect.spanLo[k] = spanLo[i];
                isect.spanHi[k] = spanHi[i];
                i++;
                k++;
            }
        }
    }

    isect.spanPtr = k;
    return isect;
}


GvsSolObjSpanList GvsSolObjSpanList::operator-() const {
    GvsSolObjSpanList slInv( spanPtr + 1 );

    GvsSurfIntersec minIntersect, maxIntersect;
    minIntersect.setDist( -FLT_MAX );
    maxIntersect.setDist(  FLT_MAX );

    if ( spanPtr == 0 ) {
        slInv.spanLo[0] = minIntersect;
        slInv.spanHi[0] = maxIntersect;
        slInv.spanPtr = 1;
    }
    else {
        int k = 0;

        if ( spanLo[0] > -FLT_MAX ) {
            slInv.spanLo[k] = minIntersect;
            slInv.spanHi[k++] = spanLo[0];
        }

        for ( int i = 1; i < spanPtr; i++ ) {
            slInv.spanLo[k] = spanHi[i-1];
            slInv.spanHi[k++] = spanLo[i];
        }

        if ( spanHi[spanPtr-1] < FLT_MAX ) {
            slInv.spanLo[k] = spanHi[spanPtr-1];
            slInv.spanHi[k++] = maxIntersect;
        }

        slInv.spanPtr = k;
    }

    return slInv;
}

