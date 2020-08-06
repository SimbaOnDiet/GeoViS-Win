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

#include "parse_transformation.h"
#include "parse_helper.h"

#include <GvsGlobalDefs.h>

#include "Parser/GvsParseScheme.h"

#include "math/TransCoordinates.h"
#include "scheme.h";

m4d::TransCoordinates *transformator;
int length;
double *oldPos;
double *newPos;
double *oldDir;
double *newDir;
m4d::vec4 oldPosP4D;
m4d::vec4 newPosP4D;
m4d::vec4 oldDirV4D;
m4d::vec4 newDirV4D;
pointer retVec, tempVec;
bool vectorTrafo = false;

//----------------------------------------------------------------------------
//         readIn
//
// List die Argumente und parst sie
//----------------------------------------------------------------------------
void readIn(scheme *sc, pointer args)
{
    vectorTrafo = false;
    if (args == sc->NIL) scheme_error("Transformation: kein Argument");

    // Initialisierung der Punkte
    transformator = new m4d::TransCoordinates();
    oldPos = new double[4];
    newPos = new double[4];

    // Bei zwei Argumenten ist es eine Vektortransformation
    if (pair_cdr(args) != sc->NIL)
    {
        vectorTrafo = true;
        oldDir = new double[4];
        newDir = new double[4];
    }

    // Liest den Punkt ein
    get_double_vec(pair_car(args), 4, oldPos, "transCartSph Punkt einlesen");

    if (vectorTrafo)
        get_double_vec(pair_car(pair_cdr(args)), 4, oldDir, "transCartSph Vektor einlesen");


    // Umwandlung in P5D-Punkt
    oldPosP4D = m4d::vec4(oldPos[0], oldPos[1], oldPos[2], oldPos[3]);

    // Umwandlung in V4D-Vektor
    if (vectorTrafo)
        oldDirV4D = m4d::vec4(oldDir[0], oldDir[1], oldDir[2], oldDir[3]);
}

//----------------------------------------------------------------------------
//         finish
//
// Umwandlung in Scheme-Format und Speicherfreigabe
//----------------------------------------------------------------------------
void finish(scheme *sc)
{
    // Umwandlung in einen Scheme-Vektor
    if (!vectorTrafo)
    {
        retVec = ((sc->vptr->mk_vector)(sc, 4));
        for (int i = 0; i < 4; i++)
        {
            tempVec = ((sc->vptr->set_vector_elem)(retVec, i, ((sc->vptr->mk_real)(sc, newPosP4D.x(i)))));
        }
    }
    else
    {
        retVec = ((sc->vptr->mk_vector)(sc, 4));
        for (int i = 0; i < 4; i++)
        {
            tempVec = ((sc->vptr->set_vector_elem)(retVec, i, ((sc->vptr->mk_real)(sc, newDirV4D.x(i)))));
        }
    }

    // Speicher freigeben
    delete oldPos;
    delete newPos;
    delete transformator;
    if (vectorTrafo)
    {
        delete oldDir;
        delete newDir;
    }
}

//----------------------------------------------------------------------------
//         gvsP_transCartSph
//----------------------------------------------------------------------------
pointer gvsP_transCartSph (scheme *sc, pointer args)
{
    readIn(sc, args);

    // Transformiert den Punkt auf Kugelkoordinaten
    if (!vectorTrafo)
        transformator->transCoordCartSph(oldPosP4D, newPosP4D);
    else
        transformator->transCoordCartSph(oldPosP4D, oldDirV4D, newPosP4D, newDirV4D);

    finish(sc);

    //Rueckgabe des transformierten Punkts an Scheme
    return retVec;
}

//----------------------------------------------------------------------------
//         gvsP_transSphCart
//----------------------------------------------------------------------------
pointer gvsP_transSphCart (scheme *sc, pointer args)
{
    readIn(sc, args);

    // Transformiert den Punkt auf Kugelkoordinaten
    if (!vectorTrafo)
        transformator->transCoordSphCart(oldPosP4D, newPosP4D);
    else
        transformator->transCoordSphCart(oldPosP4D, oldDirV4D, newPosP4D, newDirV4D);

    finish(sc);
    return retVec;
}

//----------------------------------------------------------------------------
//         gvsP_transCartCyl
//----------------------------------------------------------------------------
pointer gvsP_transCartCyl (scheme *sc, pointer args)
{
    readIn(sc, args);

    // Transformiert den Punkt auf Kugelkoordinaten
    if (!vectorTrafo)
        transformator->transCoordCartCyl(oldPosP4D, newPosP4D);
    else
        transformator->transCoordCartCyl(oldPosP4D, oldDirV4D, newPosP4D, newDirV4D);

    finish(sc);
    return retVec;
}

//----------------------------------------------------------------------------
//         gvsP_transCylCart
//----------------------------------------------------------------------------
pointer gvsP_transCylCart (scheme *sc, pointer args)
{
    readIn(sc, args);

    // Transformiert den Punkt auf Kugelkoordinaten
    if (!vectorTrafo)
        transformator->transCoordCylCart(oldPosP4D, newPosP4D);
    else
        transformator->transCoordCylCart(oldPosP4D, oldDirV4D, newPosP4D, newDirV4D);

    finish(sc);
    return retVec;
}
