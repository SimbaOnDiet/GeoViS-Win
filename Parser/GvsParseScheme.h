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
#ifndef GVS_PARSE_SCHEME_H
#define GVS_PARSE_SCHEME_H


#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <cassert>
#include <map>
#include "GvsGlobalDefs.h"
#include <Parser/parse_helper.h>


#include "scheme-private.h"
#include "scheme.h"


#ifndef T_MASKTYPE
#define T_MASKTYPE      31
#endif




enum GvsParseParamType
{
    gp_unknown = 0,
    gp_string_int,
    gp_string_double,
    gp_string_string,
    gp_string_bool,
    gp_string_gvs,
    gp_string_matrix,
    gp_string_setparamlist
};

static const std::string GvsParseParamTypeName[8] = { "unknown", "int", "double", "string", "bool", "gvs", "matrix", "splist" };

struct GvsParseParameter
{
    GvsParseParamType type;   // Parametertyp
    void* name;               // Parameter-Name: char*, symbol, ...
    int   dim;                // Dimension des Parameterwerts bei int,double
    void* val;                // Parameter-Wert: char*, int[], double[], ...
};


struct GvsParseAllowedNames
{
    GvsParseParamType  paramType;
    int                paramDim;   // fuer string verwende dim=0;
};


class GvsParseScheme
{
public:
    GvsParseScheme(scheme* s);
    GvsParseScheme(scheme* s, std::string* apNames, GvsParseAllowedNames* apType, int numNames);
    virtual ~GvsParseScheme();

    virtual pointer parse ( pointer a );  //!< Read parameters into a list

    void read_string (pointer a, const char*);

    virtual void delParamList ( void );   //!< Clear parameter list


    int   getNumParam  ( void );
    bool  getParameter ( const char* name, char*   sval );
    bool  getParameter ( const char* name, std::string  &sval );
    bool  getParameter ( const char* name, std::string  &sval, int &num);
    bool  getParameter ( const char* name, int     &ival);
    bool  getParameter ( const char* name, int*    ival );
    bool  getParameter ( const char* name, double  &dval);
    bool  getParameter ( const char* name, double* dval );
    bool  getParameter ( const char* name, double* dval, int &num, bool fixpos = true );

    bool  getParameter ( const char* name, bool    &bval );

    bool  getParameter ( const char* name, m4d::Matrix<double,2,3>* mat );
    bool  getParameter ( const char* name, m4d::Matrix<double,3,4>* mat );

    bool  getParameter ( const char* name, GvsSetParamList* splist, int &num );

    std::string        getParamName ( const int k ) const;
    int                getParamDim  ( const int k ) const;
    int                getParamDim  ( const char* name ) const;
    GvsParseParamType  getParamType ( const int k ) const;
    GvsParseParamType  getParamType ( const char* name ) const;


    // Suche in Parameterliste nach dem Parameter 'name'
    bool  testArgType  ( const char* name, int &num );

    void               setAllowedName      ( const std::string name, const GvsParseParamType pType, const int dim);
    GvsParseParamType  getAllowedParamType ( const std::string name );
    int                getAllowedParamDim  ( const std::string name );

    bool  allowedName       ( const char*  name );
    void  testParamNames    ( const char* name );
    void  printAllowedNames ( );

    bool  transType ( const GvsDataType datType, GvsParseParamType &pType, int &dim );

    void  print () const;

private:
    scheme* sc;      //!< pointer to scheme
    pointer args;    //!< pointer to current argument

    int  numParameter;  //!< number of parameters
    std::vector<GvsParseParameter*> paramListe;   //!< parameter list

    std::map<std::string,GvsParseAllowedNames>           allowedParamNames;
    std::map<std::string,GvsParseAllowedNames>::iterator allowedParamNamesPtr;
};


inline int GvsParseScheme :: getNumParam ( void ) {
    return numParameter;
}


#endif
