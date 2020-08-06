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
/**
    Translation
    @verbatim
    (translate #(double double double))@endverbatim

    Transformation matrices can also be composited. For example, two rotations could read like this:
    @verbatim
    (rotate-obj "z" 20.0 (rotate-obj "x" 10.0))@endverbatim
    This can, in principle, be repeated indefinitely often.
*/

#include "parse_transformObject.h"
#include "parse_helper.h"
#include "Parser/GvsParseScheme.h"

#include <GvsGlobalDefs.h>
#include "scheme.h"

#include "m4dGlobalDefs.h"
#include "math/TransfMat.h"

int dim = 0;
double *transformVec;
m4d::vec2 dir2d;
m4d::vec3 dir3d;
std::string achse;

m4d::Matrix<double,2,3> *inMatrix2d;
m4d::Matrix<double,3,4> *inMatrix3d;
m4d::Matrix<double,2,3> *transformMat2d;
m4d::Matrix<double,3,4> *transformMat3d;

pointer matPointer;

/**
 * @brief gvsP_translateObj
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_translateObj (scheme *sc, pointer args)
{
    if (args == sc->NIL) scheme_error("Transformation translate-object: kein Argument");

    // Bestimme als erstes die Dimension der Transformation mit Hilfe des
    // Translationsvektors
    dim = (sc->vptr->vector_length)( pair_car(args) );
    if ( !((dim ==2) || (dim == 3)) )
        scheme_error("Transformation translate-object: Falsche Dimension");

    transformVec = new double[dim];

    // Parsen der Argumente
    get_double_vec(pair_car(args), dim, transformVec, "Transformation translate-object: Vektor einlesen");

    matPointer = sc->NIL;
    if ( dim == 2 )
    {
        // Erzeuge die Transformationsmatrix
        transformMat2d = new m4d::TranslateMat2D(transformVec[0], transformVec[1]);

        // Weiteres Parsen falls verschachtelte Transformation
        if (pair_cdr(args) != sc->NIL)
        {
            inMatrix2d = new m4d::Matrix<double,2,3>();
            get_matrix(pair_car(pair_cdr(args)), inMatrix2d, "Transformation translate-object: Matrix einlesen");

            // Multipliziere sie mit der eingelesenen Matrix
            //*inMatrix = *(dynamic_cast<m4d::Matrix<double,3,4>*>(translateMat3d)) * *inMatrix;
            *inMatrix2d = *transformMat2d * *inMatrix2d;
            // Umwandlung
            mk_matrix( matPointer, inMatrix2d, "Transformation translate-object: Matrix fuer Scheme bauen");
            delete inMatrix2d;
        }
        else
        {
            // Umwandlung
            mk_matrix( matPointer, transformMat2d, "Transformation translate-object: Matrix fuer Scheme bauen");
        }
        //delete transformMat2d;
    }
    else
    {
        // Erzeuge die Transformationsmatrix
        transformMat3d = new m4d::TranslateMat3D(transformVec[0], transformVec[1], transformVec[2]);

        // Weiteres Parsen falls verschachtelte Transformation
        if (pair_cdr(args) != sc->NIL)
        {
            inMatrix3d = new m4d::Matrix<double,3,4>();
            get_matrix(pair_car(pair_cdr(args)), inMatrix3d, "Transformation translate-object: Matrix einlesen");

            // Multipliziere sie mit der eingelesenen Matrix
            //*inMatrix = *(dynamic_cast<m4d::Matrix<double,3,4>*>(translateMat3d)) * *inMatrix;
            *inMatrix3d = *transformMat3d * *inMatrix3d;
            // Umwandlung
            mk_matrix( matPointer, inMatrix3d, "Transformation translate-object: Matrix fuer Scheme bauen");
            delete inMatrix3d;
        }
        else
        {
            // Umwandlung
            mk_matrix( matPointer, transformMat3d, "Transformation translate-object: Matrix fuer Scheme bauen");
        }
        delete transformMat3d;
    }

    delete transformVec;

    return matPointer;
}

/**
 * Constructs a rotation matrix from an axis label or an axis vector together with an angle in degrees.
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_rotateObj (scheme *sc, pointer args) {
    double angle;
    //std::cerr << "rotateObj.....................\n";
    if (args == sc->NIL) scheme_error("Transformation rotate-object: no argument");

    // If there is an axis label, then it is a three-dimensional rotation.
    if (is_string(pair_car(args))) {
        dim = 3;
        transformVec = new double[dim];
        get_string(pair_car(args), achse, "Transformation rotate-obj: axis label should be a string!");

        if ( (achse == "x") || (achse == "X") )      dir3d = m4d::vec3(1.0, 0.0, 0.0);
        else if ( (achse == "y") || (achse == "Y") ) dir3d = m4d::vec3(0.0, 1.0, 0.0);
        else if ( (achse == "z") || (achse == "Z") ) dir3d = m4d::vec3(0.0, 0.0, 1.0);
        else scheme_error("Transformation rotate-object: wrong axis label");

        get_double(pair_car(pair_cdr(args)), &angle, "Transformation rotate-object: read angle");
    }

    // Ist das erste Argument ein Vektor, so muss die Dimension bestimmt werden.
    // Es kann sich in diesem Fall sowohl um eine zwei- als auch um eine drei-
    // dimensionale Transformation handeln.
    else if ( is_vector(pair_car(args)) )
    {
        dim = (sc->vptr->vector_length)( pair_car(args) );
        if ( !((dim ==2) || (dim == 3)) )
            scheme_error("Transformation rotate-object: wrong dimension");

        transformVec = new double[dim];
        get_double_vec(pair_car(args), dim, transformVec, "Transformation rotate-object: read vector");
        if ( dim == 3 )
            dir3d = m4d::vec3(transformVec[0], transformVec[1], transformVec[2]);

        get_double(pair_car(pair_cdr(args)), &angle, "Transformation rotate-object: read angle");
    }

    // Ist das erste Argument eine Zahl, also der Rotationswinkel, so muss die
    // Transformation zweidimensional sein.
    else if ( (is_real(pair_car(args))) || (is_integer(pair_car(args))) )
    {
        dim = 2;
        transformVec = new double[dim];
        transformVec[0] = 0.0;
        transformVec[1] = 0.0;
        // Rotationswinkel
        get_double(pair_car(args), &angle, "Transformation rotate-object: Winkel einlesen");
    }

    if ( dim == 2)
    {
        // Construct 2D transformation matrix
        transformMat2d = new m4d::RotateMat2D(transformVec[0], transformVec[1], angle * GVS_deg2rad);

        // Bei drei Argumenten ist es eine geschachtelte Transformation
        if (pair_cdr(pair_cdr(args)) != sc->NIL)
        {
            inMatrix2d = new m4d::Matrix<double,2,3>();

            get_matrix(pair_car(pair_cdr(pair_cdr(args))), inMatrix2d, "Transformation rotate-object: Matrix einlesen");

            // Multipliziere sie mit der eingelesenen Matrix
            //*inMatrix = *(dynamic_cast<m4d::Matrix<double,3,4>*>(transformMat)) * *inMatrix;
            *inMatrix2d = *transformMat2d * *inMatrix2d;

            // Umwandlung
            mk_matrix( matPointer, inMatrix2d, "Transformation rotate-object: Matrix fuer Scheme bauen");
            delete inMatrix2d;
        }
        else
        {
            // Umwandlung
            mk_matrix( matPointer, transformMat2d, "Transformation rotate-object: Matrix fuer Scheme bauen");
        }
        //  delete transformMat2d;
    }
    else {
        // Construct 3D transformation matrix
        transformMat3d = new m4d::RotateMat3D(dir3d, angle * GVS_deg2rad);
        //transformMat3d->print();

        // Bei drei Argumenten ist es eine geschachtelte Transformation
        if (pair_cdr(pair_cdr(args)) != sc->NIL)
        {
            inMatrix3d = new m4d::Matrix<double,3,4>();

            get_matrix(pair_car(pair_cdr(pair_cdr(args))), inMatrix3d, "Transformation rotate-object: Matrix einlesen");

            // Multipliziere sie mit der eingelesenen Matrix
            //*inMatrix = *(dynamic_cast<m4d::Matrix<double,3,4>*>(transformMat)) * *inMatrix;
            *inMatrix3d = *transformMat3d * *inMatrix3d;

            mk_matrix( matPointer, inMatrix3d, "Transformation rotate-object: Matrix fuer Scheme bauen");
            delete inMatrix3d;
        }
        else {
            mk_matrix( matPointer, transformMat3d, "Transformation rotate-object: Matrix fuer Scheme bauen");
            // std::cerr << "length: " << (matPointer->_flag) << std::endl;
        }
        //   delete transformMat3d;
    }

    delete transformVec;
    return matPointer;
}

/**
 * @brief gvsP_scaleObj
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_scaleObj (scheme *sc, pointer args) {
    if (args == sc->NIL) {
        scheme_error("Transformation scale-obj: no argument");
    }

    // Bestimme als erstes die Dimension der Transformation ueber die
    // Dimension des uebergebenen Vektors
    dim = (sc->vptr->vector_length)( pair_car(args) );
    if ( !((dim ==2) || (dim == 3)) ) {
        scheme_error("Transformation scale-obj: wrong dimension");
    }

    transformVec = new double[dim];
    get_double_vec(pair_car(args), dim, transformVec, "Transformation scale-obj: read vector");

    if ( dim == 2 ) {
        // generate transformation matrix
        transformMat2d = new m4d::ScaleMat2D(transformVec[0], transformVec[1]);

        // Weiteres Parsen falls verschachtelte Transformation
        if (pair_cdr(args) != sc->NIL)
        {
            inMatrix2d = new m4d::Matrix<double,2,3>();
            get_matrix(pair_car(pair_cdr(args)), inMatrix2d, "Transformation scale-obj: Matrix einlesen");

            // Multipliziere sie mit der eingelesenen Matrix
            //*inMatrix = *(dynamic_cast<m4d::Matrix<double,3,4>*>(translateMat3d)) * *inMatrix;
            *inMatrix2d = *transformMat2d * *inMatrix2d;
            // Umwandlung
            mk_matrix( matPointer, inMatrix2d, "Transformation scale-obj: Matrix fuer Scheme bauen");
            delete inMatrix2d;
        }
        else
        {
            // Umwandlung
            mk_matrix( matPointer, transformMat2d, "Transformation scale-obj: Matrix fuer Scheme bauen");
        }
        delete transformMat2d;
    }
    else
    {
        // Erzeuge die Transformationsmatrix
        transformMat3d = new m4d::ScaleMat3D(transformVec[0], transformVec[1], transformVec[2]);

        // Weiteres Parsen falls verschachtelte Transformation
        if (pair_cdr(args) != sc->NIL)
        {
            inMatrix3d = new m4d::Matrix<double,3,4>();
            get_matrix(pair_car(pair_cdr(args)), inMatrix3d, "Transformation scale-obj: Matrix einlesen");

            // Multipliziere sie mit der eingelesenen Matrix
            //*inMatrix = *(dynamic_cast<m4d::Matrix<double,3,4>*>(translateMat3d)) * *inMatrix;
            *inMatrix3d = *transformMat3d * *inMatrix3d;
            // Umwandlung
            mk_matrix( matPointer, inMatrix3d, "Transformation scale-obj: Matrix fuer Scheme bauen");
            delete inMatrix3d;
        }
        else
        {
            // Umwandlung
            mk_matrix( matPointer, transformMat3d, "Transformation scale-obj: Matrix fuer Scheme bauen");
        }
        //delete transformMat3d;
    }

    delete [] transformVec;

    return matPointer;
}
