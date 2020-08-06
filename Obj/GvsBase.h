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
//DONE READING
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

#ifndef GVS_BASE_H
#define GVS_BASE_H

//
//  A parameter can be made available as ChangeParameter, if there is a corresponding
//  SetParam implementation in the derived class and an 'AddParam' was set in the
//  constructor
//

#include <iostream>
#include <algorithm>
#include <map>

#include "GvsGlobalDefs.h"
#include "metric/m4dMetric.h"

#include "m4dGlobalDefs.h"
class GvsSceneObj;

/**
 * @brief The GvsBase class
 */
class GvsBase
{
public:
    GvsBase();
    virtual ~GvsBase();

    virtual bool  AddParam ( std::string pName, const GvsDataType type );
    virtual void  DelParam ( std::string pName );
    virtual void  DelAllParam( );

    virtual bool  SetParam ( std::string pName, int val );
    virtual bool  GetParam ( std::string pName, int &val );

    virtual bool  SetParam ( std::string pName, double val );
    virtual bool  GetParam ( std::string pName, double &val );

    virtual bool  SetParam ( std::string pName, m4d::ivec2 vec );
    virtual bool  GetParam ( std::string pName, m4d::ivec2 &vec );

    virtual bool  SetParam ( std::string pName, m4d::ivec3 vec );
    virtual bool  GetParam ( std::string pName, m4d::ivec3 &vec );

    virtual bool  SetParam ( std::string pName, m4d::ivec4 vec );
    virtual bool  GetParam ( std::string pName, m4d::ivec4 &vec );

    virtual bool  SetParam ( std::string pName, m4d::vec2 pt );
    virtual bool  GetParam ( std::string pName, m4d::vec2 &pt );

    virtual bool  SetParam ( std::string pName, m4d::vec3 pt );
    virtual bool  GetParam ( std::string pName, m4d::vec3 &pt );

    virtual bool  SetParam ( std::string pName, m4d::vec4 pt );
    virtual bool  GetParam ( std::string pName, m4d::vec4 &pt );

    virtual bool  SetParam ( std::string pName, m4d::Matrix<double,2,3> mat );
    virtual bool  GetParam ( std::string pName, m4d::Matrix<double,2,3> &mat );

    virtual bool  SetParam ( std::string pName, m4d::Matrix<double,3,4> mat );
    virtual bool  GetParam ( std::string pName, m4d::Matrix<double,3,4> &mat );

    virtual bool  SetParam ( std::string pName, std::string txt );
    virtual bool  GetParam ( std::string pName, std::string &txt );


    bool  IsValidParamName ( std::string pName );
    bool  IsValidParam ( std::string pName, GvsDataType dataType );


    virtual unsigned int GetNumParams () const;
    virtual std::string  GetParamName ( unsigned int i ) const;
    virtual GvsDataType  GetDataType  ( std::string pName );
    virtual void         PrintAllParameter();


    virtual void  Add ( GvsSceneObj* obj );
    virtual void  Print ( FILE* fptr = stderr ) = 0;

    void        lowCase    ( std::string &s );
    std::string getLowCase ( std::string &s ) const;

    void ResetID() { mID = 0; }
    int  GetID() { return mID; }

protected:
    unsigned int mNumParam;
    std::map<std::string,gvs_parameter>            mParam;
    std::map<std::string,gvs_parameter>::iterator  mParamPtr;

    int mID;
    static int mObjCounter;
};


/**
 * @brief The Gvsm4dMetricDummy class
 *
 *  Necessary to subclass a metric from GvsBase !!
 */
class Gvsm4dMetricDummy : public GvsBase
{
public:
    Gvsm4dMetricDummy (m4d::Metric* cMetric);
    virtual ~Gvsm4dMetricDummy();
    m4d::Metric* m4dMetric;
    virtual bool SetParam ( std::string pName, double val);
    virtual bool GetParam ( std::string pName, double &val);
    virtual void Print( FILE* fptr = stderr );
};

#endif
