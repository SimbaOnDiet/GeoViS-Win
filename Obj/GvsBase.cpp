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

#include "GvsBase.h"

int GvsBase::mObjCounter = 0;

GvsBase::GvsBase() :
    mNumParam(0) {
}

GvsBase::~GvsBase() {
    if (!mParam.empty()) {
        mParam.clear();
    }
    mNumParam = 0;
}

bool GvsBase::AddParam( std::string pName, const GvsDataType type ) {
    // parameter names are stored only in lower case format
    lowCase(pName);

    // test, if parameter already exists
    mParamPtr = mParam.find(pName);
    if (mParamPtr == mParam.end()) {
        gvs_parameter par = {type,0};
        mParam.insert(std::pair<std::string,gvs_parameter>(pName,par));
        //mParam.insert(std::pair<std::string,gvs_parameter>(pName,(gvs_parameter){type,0}));
        mNumParam++;
        return true;
    } else {
        fprintf(stderr,"Parameter %s already exists!\n",pName.c_str());
        return false;
    }
}

void GvsBase::DelParam( std::string pName ) {
    lowCase(pName);

    mParamPtr = mParam.find(pName);
    if (mParamPtr == mParam.end()) {
        fprintf(stderr,"Parameter %s not available, hence cannot be deleted!\n",pName.c_str());
    } else {
        gvs_parameter par = mParamPtr->second;
        switch (par.type) {
            default:
                break;
            case gvsDT_INT:
                delete (int*)par.val;
                break;
            case gvsDT_DOUBLE:
                delete (double*)par.val;
                break;
            case gvsDT_VEC3:
                delete (m4d::vec3*)par.val;
                break;
            case gvsDT_VEC4:
                delete (m4d::vec4*)par.val;
                break;
            case gvsDT_STRING:
                delete (std::string*)par.val;
                break;
        }

        mParam.erase(mParamPtr);
        mNumParam = mParam.size();
    }
}

void GvsBase::DelAllParam() {
    if (!mParam.empty()) {

        std::vector<std::string> plist;
        mParamPtr = mParam.begin();
        do {
            plist.push_back(mParamPtr->first);
            mParamPtr++;
        } while (mParamPtr!=mParam.end());

        for(unsigned int i=0; i<plist.size(); i++) {
            DelParam(plist[i]);
        }
        mParam.clear();
    }
    mNumParam = 0;
}


bool GvsBase::SetParam( std::string pName, int val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete static_cast<int*>((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_INT) {
        fprintf(stderr,"Parameter %s expects value of type %s !\n",
                pName.c_str(),GvsDataTypeName[(mParamPtr->second).type].c_str());
        return false;
    }
    int* value = new int(val);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam( std::string pName, double val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (double*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_DOUBLE) {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    double* value = new double(val);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam( std::string pName, m4d::vec2 pt ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::vec2*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_VEC2) {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    m4d::vec2* value = new m4d::vec2(pt);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam ( std::string pName, m4d::ivec2 vec ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::ivec2*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_IVEC2) {
        fprintf(stderr,"SetParam(): %s expected %s.\n",pName.c_str(),GvsDataTypeName[mParamPtr->second.type].c_str());
        return false;
    }
    m4d::ivec2* value = new m4d::ivec2(vec);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam ( std::string pName, m4d::ivec3 vec ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::ivec3*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_IVEC3) {
        fprintf(stderr,"SetParam(): %s expected %s.\n",pName.c_str(),GvsDataTypeName[mParamPtr->second.type].c_str());
        return false;
    }
    m4d::ivec3* value = new m4d::ivec3(vec);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam ( std::string pName, m4d::ivec4 vec ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::ivec4*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_IVEC4) {
        fprintf(stderr,"SetParam(): %s expected %s.\n",pName.c_str(),GvsDataTypeName[mParamPtr->second.type].c_str());
        return false;
    }
    m4d::ivec4* value = new m4d::ivec4(vec);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam( std::string pName, m4d::vec3 pt ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::vec3*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_VEC3) {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    m4d::vec3* value = new m4d::vec3(pt);
    (mParamPtr->second).val = value;
    return true;
}

bool GvsBase::SetParam( std::string pName, m4d::vec4 pt ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::vec4*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_VEC4) {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    m4d::vec4* value = new m4d::vec4(pt);
    (mParamPtr->second).val = value;
    return true;
}


bool GvsBase::SetParam ( std::string pName, m4d::Matrix<double,2,3> mat ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::Matrix<double,2,3>*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_MAT2D) {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    m4d::Matrix<double,2,3>* value = new m4d::Matrix<double,2,3>(mat);
    (mParamPtr->second).val = value;
    return true;
}


bool GvsBase::SetParam ( std::string pName, m4d::Matrix<double,3,4> mat ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (m4d::Matrix<double,3,4>*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_MAT3D)
    {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    m4d::Matrix<double,3,4>* value = new m4d::Matrix<double,3,4>(mat);
    (mParamPtr->second).val = value;
    return true;
}


bool GvsBase::SetParam( std::string pName, std::string txt ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    delete (std::string*)((mParamPtr->second).val);

    if ((mParamPtr->second).type != gvsDT_STRING) {
        std::cerr << "setParam(): " << pName << " erwartet " << GvsDataTypeName[(mParamPtr->second).type] << std::endl;
        return false;
    }
    std::string* value = new std::string(txt);
    (mParamPtr->second).val = value;
    return true;
}


bool GvsBase::GetParam( std::string pName, int &val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_INT) {
        std::cerr << "getParam() ... integer erwartet!" << std::endl;
        return false;
    }
    val = *(static_cast<int*>((mParamPtr->second).val));
    return true;
}

bool GvsBase::GetParam( std::string pName, double &val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_DOUBLE) {
        std::cerr << "getParam() ... double erwartet!" << std::endl;
        return false;
    }
    val = *(static_cast<double*>((mParamPtr->second).val));
    return true;
}


bool GvsBase::GetParam( std::string pName, m4d::vec3 &val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_VEC3) {
        std::cerr << "setParam() ... m4d::vec2 erwartet!" << std::endl;
        return false;
    }
    val = *((m4d::vec3*)(mParamPtr->second).val);
    return true;
}

bool GvsBase::GetParam( std::string pName, m4d::vec2 &val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_VEC2) {
        std::cerr << "setParam() ... m4d::vec2 erwartet!" << std::endl;
        return false;
    }
    val = *((m4d::vec2*)(mParamPtr->second).val);
    return true;
}

bool GvsBase::GetParam ( std::string pName, m4d::ivec2 &vec ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    if ((mParamPtr->second).type != gvsDT_IVEC2) {
        fprintf(stderr,"GetParam():  m4d::ivec2 expected!\n");
        return false;
    }
    vec = *(static_cast<m4d::ivec2*>((mParamPtr->second).val));
    return true;
}

bool GvsBase::GetParam ( std::string pName, m4d::ivec3 &vec ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    if ((mParamPtr->second).type != gvsDT_IVEC3) {
        fprintf(stderr,"GetParam():  m4d::ivec3 expected!\n");
        return false;
    }
    vec = *(static_cast<m4d::ivec3*>((mParamPtr->second).val));
    return true;
}

bool GvsBase::GetParam ( std::string pName, m4d::ivec4 &vec ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }
    if ((mParamPtr->second).type != gvsDT_IVEC4) {
        fprintf(stderr,"GetParam():  m4d::ivec4 expected!\n");
        return false;
    }
    vec = *(static_cast<m4d::ivec4*>((mParamPtr->second).val));
    return true;
}

bool GvsBase::GetParam( std::string pName, m4d::vec4 &val ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_VEC4) {
        std::cerr << "setParam() ... m4d::vec4 erwartet!" << std::endl;
        return false;
    }
    val = *((m4d::vec4*)(mParamPtr->second).val);
    return true;
}

bool GvsBase::GetParam ( std::string pName, m4d::Matrix<double,2,3> &mat ) {
    lowCase(pName);
    mParamPtr = mParam.find(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_MAT2D) {
        std::cerr << "setParam() ... m4d::vec3 erwartet!" << std::endl;
        return false;
    }

    mat = *((m4d::Matrix<double,2,3>*)(mParamPtr->second).val);
    return true;
}

//----------------------------------------------------------------------------
//         getParam  TMat3D
//----------------------------------------------------------------------------
bool GvsBase::GetParam ( std::string pName, m4d::Matrix<double,3,4> &mat ) {
    lowCase(pName);
    mParamPtr = mParam.find(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_MAT3D) {
        std::cerr << "setParam() ... m4d::vec3 erwartet!" << std::endl;
        return false;
    }

    mat = *((m4d::Matrix<double,3,4>*)(mParamPtr->second).val);
    return true;
}


bool GvsBase::GetParam( std::string pName, std::string &txt ) {
    lowCase(pName);
    if (!IsValidParamName(pName)) {
        return false;
    }

    if ((mParamPtr->second).type != gvsDT_STRING) {
        std::cerr << "setParam() ... STRING erwartet!" << std::endl;
        return false;
    }
    txt = *((std::string*)(mParamPtr->second).val);
    return true;
}


bool GvsBase::IsValidParamName( std::string pName ) {
    mParamPtr = mParam.find(pName);
    if (mParamPtr == mParam.end()) {
        fprintf(stderr,"Parameter %s not available\n",pName.c_str());
        return false;
    }
    return true;
}

bool GvsBase::IsValidParam ( std::string pName, GvsDataType dataType ) {
    if (!IsValidParamName(pName)) {
        return false;
    }
    if ((mParamPtr->second).type != dataType) {
        std::cerr << "setParam() ... " << GvsDataTypeName[(mParamPtr->second).type] << " erwartet!" << std::endl;
        return false;
    }
    return true;
}

unsigned int GvsBase::GetNumParams() const {
    return mNumParam;
}

std::string GvsBase::GetParamName( unsigned int i ) const {
    if (i>=mNumParam) {
        return std::string();
    }
    std::string pName;

    unsigned int k=0;
    std::map<std::string,gvs_parameter>          par = mParam;
    std::map<std::string,gvs_parameter>::iterator parPtr = par.begin();
    do {
        pName = (parPtr->first);
        parPtr++;
        k++;
    } while ((k<=i) && (parPtr!=par.end()));
    return pName;
}

GvsDataType GvsBase::GetDataType( std::string pName ) {
    if (!IsValidParamName(pName)) {
        return gvsDT_UNKNOWN;
    }
    return (mParamPtr->second).type;
}

void GvsBase::PrintAllParameter() {
    mParamPtr = mParam.begin();
    do {
        std::cerr << (mParamPtr->first) << std::endl;
        mParamPtr++;
    } while (mParamPtr!=mParam.end());
}


void GvsBase::Add ( GvsSceneObj*  ) {
    std::cerr << "Gvs :: add()  not implemented here!" << std::endl;
}


//----------------------------------------------------------------------------
//       lowCase
//----------------------------------------------------------------------------
void
GvsBase :: lowCase(std::string &s)
{
  int l = s.length();
  for (int i=0;i<l;i++) s[i] = tolower(s[i]);
  // transform(s.begin(),s.end(),s.begin(),tolower);
}

std::string
GvsBase :: getLowCase ( std::string &s ) const
{
  std::string sOut = s;
  int l = sOut.length();
  for (int i=0;i<l;i++)
  {
    tolower(sOut[i]);
  }

  return sOut;
}



Gvsm4dMetricDummy::Gvsm4dMetricDummy (m4d::Metric* cMetric) {
    m4dMetric=cMetric;
}

bool Gvsm4dMetricDummy::SetParam( std::string pName, double val ) {
    return m4dMetric->setParam(pName, val);
}


bool Gvsm4dMetricDummy::GetParam( std::string pName, double &val ) {
    return m4dMetric->getParam(pName, val);
}

void Gvsm4dMetricDummy::Print( FILE* fptr ) {
    m4dMetric->print(fptr);
    fprintf(fptr,"\n");
}

Gvsm4dMetricDummy::~Gvsm4dMetricDummy () {
   //delete m4dMetric;
   //GvsBase::~GvsBase ();
}
