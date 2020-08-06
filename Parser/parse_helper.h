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
#ifndef GVS_PARSE_HELPER_H
#define GVS_PARSE_HELPER_H

extern "C" {
#include "scheme-private.h"
#include "scheme.h"
}


#include "Obj/GvsBase.h"
#include "Shader/GvsShader.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "metric/m4dMetric.h"

class GvsParseScheme;

extern scheme sc;

enum GvsType
{
    gtunknown = 0,
    gtCamera,
    gtDevice,
    gtGeodSolver,
    gtLightMgr,
    gtLight,
    gtLocTed,
    gtMetric,
    gtMotion,
    gtProjector,
    gtRayGen,
    gtSceneObj,
    gtShader,
    gtTexture
};

const std::string GvsTypeName[14] = {
    "gtUnknown",
    "gtCamera",
    "gtDevice",
    "gtGeodSolver",
    "gtLightMgr",
    "gtLight",
    "gtLocTed",
    "gtMetric",
    "gtMotion",
    "gtProjector",
    "gtRayGen",
    "gtSceneObj",
    "gtShader",
    "gtTexture"
};

struct GvsTypeID
{
    //string  name;
    GvsType  gvsType;
    int      vectorID;
    GvsBase* gvsObject;
};


struct GvsSetParamList    // setparameter-Liste
{
    char*          objectIDname;
    char*          paramName;
    gvs_parameter  param;
};



void scheme_error ( const std::string& msg );

// ----------- String-Funktionen -------------
void        lowCase    ( std::string &s );
std::string getLowCase ( std::string s );
void        appendNum  ( std::string &s, int num );

int searchCharArray   ( char**  stringArray, const int num, const char* searchItem );
int searchStringArray ( const std::string* stringArray, const int num, const char* searchItem );

std::string getType ( pointer arg );  //!< type of scheme element

pointer gvsP_print ( scheme *sc, pointer args );
pointer gvsP_m4d_metriclist ( scheme *sc, pointer args );
pointer gvsP_m4d_solverlist ( scheme *sc, pointer args );
pointer gvsP_set_parameter ( scheme *sc, pointer args );
pointer gvsP_get_parameter ( scheme *sc, pointer args );
pointer gvsP_exit( scheme *sc, pointer args );
pointer gvsP_getenv( scheme *sc, pointer args );

std::map<std::string,GvsTypeID>::iterator getIDptr ( GvsParseScheme* gpSc,
                                                     const std::string objName,
                                                     const std::string suchName,
                                                     const char* token,
                                                     const GvsType gtype,
                                                     int num = 0);

std::map<std::string,GvsTypeID>::iterator getOneOfIDptr ( GvsParseScheme* gpSc,
                                                          const std::string objName,
                                                          const std::string suchName,
                                                          const char* token,
                                                          const GvsType gtype1, const GvsType gtype2, GvsType &gtypeFound,
                                                          int num = 0);


m4d::Metric* readMetric (const std::string name, GvsParseScheme* gP );

GvsShader* readShader ( const std::string name, GvsParseScheme* gP );
void  readMetric ( const std::string name, GvsParseScheme* gP, m4d::Metric* currMetric );

void get_double ( pointer s_x, double* x, const std::string msg = "" );
void get_int    ( pointer s_x, int* x,    const std::string msg = "" );
void get_char   ( pointer s_c, char* tx,  const std::string msg = "" );
void get_string ( pointer s_s, std::string& s, const std::string msg = "" );

void get_int_vec    ( pointer s_vec, int dim, int p[], const std::string msg = "" );
void get_double_vec ( pointer s_vec, int dim, double p[], const std::string msg = "" );

void get_matrix_size ( pointer s_vec, int &n, int &m );

// void get_matrix(pointer s_vec, Matrix<mType, n, m> *mat, const string msg = "");


/**
 *   Read a transformation matrix from scheme
 */
template <class mType, int n, int m> void get_matrix(pointer s_vec, m4d::Matrix<mType, n, m> *mat, const std::string msg)
{
    double *tempV = new double[m];
    if (!is_vector(s_vec)) scheme_error(msg + ": get_transformation_matrix: Falsches Format der Matrix");

    for (int i = 0; i < n; i++)
    {
        get_double_vec( ((*sc.vptr->vector_elem)(s_vec, i)) , m, tempV, "get_transformation_matrix: Matrix einlesen");
        for (int j = 0; j < m; j++)
        {
            mat->setCoeff(i, j, tempV[j]);
        }
    }
    delete [] tempV;
}


/**
 * Convert a m4d matrix into a scheme pointer
 */
template <class mType, int n, int m> void mk_matrix(pointer &retMat, m4d::Matrix<mType, n, m> *inMat, const std::string ) {
    // column vector of row vectors...
    pointer vecRowPointer;
    UNUSED_ATTR pointer tempPointer;

    retMat = ((*sc.vptr->mk_vector)( &sc, n ));
    for (int i = 0; i < n; i++) {
        vecRowPointer = ((*sc.vptr->mk_vector)(&sc, m));
        for (int j = 0; j < m; j++) {
            tempPointer = ((*sc.vptr->set_vector_elem)(vecRowPointer, j, ((*sc.vptr->mk_real)( &sc, inMat->getCoeff(i, j) ))));
        }
        tempPointer = ((*sc.vptr->set_vector_elem)(retMat, i, vecRowPointer));
    }
}

pointer gvsP_vec3 (scheme *sc, pointer args);
pointer gvsP_vec4 (scheme *sc, pointer args);

#endif
