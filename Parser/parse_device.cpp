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
        (init-device '(type "standard")
                     '(camera "camID")
                     '(projector "projID")
                     '(metric "metricID")
                     '(light-mgr "lmID")
                    [ '(setparam ("metricID" "mass" 0.3)) ]
        )
 */

#include "parse_device.h"
#include "parse_helper.h"
#include "GvsParseScheme.h"

#include "Dev/GvsDevice.h"
#include "Light/GvsLightSrcMgr.h"
#include "Utils/GvsGeodSolver.h"

#include "scheme-private.h"

extern std::vector<GvsCamera*>          gpCamera;
extern std::vector<GvsProjector*>       gpProjector;

extern std::vector<Gvsm4dMetricDummy*>  gpMetric;
extern std::vector<GvsGeodSolver*>      gpSolver;
extern std::vector<GvsSceneObj*>        gpSceneObj;
extern std::vector<GvsLightSrcMgr*>     gpLightMgr;
extern std::vector<GvsDevice*>          gpDevice;

extern std::map<std::string,GvsTypeID>  gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator  gpTypeIDptr;


//----------------------------------------------------------------------------
//         gvsP_init_device
//----------------------------------------------------------------------------
pointer gvsP_init_device (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_init_device..........\n";
#endif
    if (args == sc->NIL) scheme_error("init-device: no arguments");
    if (!is_pair(args)) scheme_error("init-device: less arguments");

    std::string allowedNames[] = {"type","camera","projector","metric","lightmgr","obj","id","setparam"};
    GvsParseAllowedNames allowedTypes[] = {{gp_string_string,0}, // type
                                           {gp_string_string,0}, // camera
                                           {gp_string_string,0}, // projector
                                           {gp_string_string,0}, // metric
                                           {gp_string_string,0}, // light-mgr
                                           {gp_string_string,0}, // obj
                                           {gp_string_string,0}, // id
                                           {gp_string_setparamlist,0} // setparam
                                          };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,8);
    args = gvsParser->parse(args);
    gvsParser->testParamNames("init-device");

    GvsDevice* currDevice = new GvsDevice();

    std::string msg;

    std::string deviceName;
    gvsParser->getParameter("type",deviceName);

    if (deviceName=="standard" || deviceName=="standard_manual") {

        if (!gpMetric.empty()) {
            currDevice->metric = gpMetric[0]->m4dMetric;
        }
        else {
            scheme_error("init-device: metric is missing!");
        }

        if (!gpCamera.empty()) {
            currDevice->camera  = gpCamera[0];
        }
        else {
            scheme_error("init-device: camera is missing!");
        }

        if (!gpProjector.empty()) {
            currDevice->projector = gpProjector[0];
        }
        else {
            scheme_error("init-device: projector is missing!");
        }

        if (!gpLightMgr.empty()) {
            currDevice->lightSrcMgr = gpLightMgr[0];
        }
        else {
            scheme_error("init-device: light source manager is missing!");
        }

        if (!gpSceneObj.empty()) {
            currDevice->sceneGraph = gpSceneObj[gpSceneObj.size()-1];
        }
        else {
            scheme_error("init-device: scene object is missing!");
        }
    }

    if (deviceName=="standard_manual") {
        currDevice->setManual(true);
    }

    //
    // -------- set specific camera
    //
    if (gpCamera.size()>1)
    {
        gpTypeIDptr = getIDptr(gvsParser,"Device","Camera","camera",gtCamera);
        currDevice->camera = gpCamera[(gpTypeIDptr->second).vectorID];
    }

    //
    // -------- set specific projector
    //
    if (gpProjector.size()>1)
    {
        gpTypeIDptr = getIDptr(gvsParser,"Device","Projector","projector",gtProjector);
        currDevice->projector = gpProjector[(gpTypeIDptr->second).vectorID];
    }

    //
    // -------- set specific metric
    //

    if (gpMetric.size()>1) {
        gpTypeIDptr = getIDptr(gvsParser,"Device","Metric","metric",gtMetric);
        currDevice->metric = gpMetric[(gpTypeIDptr->second).vectorID]->m4dMetric;
    }

    //
    // -------- set specific light source manager
    //
    if (gpLightMgr.size()>1)
    {
        gpTypeIDptr = getIDptr(gvsParser,"Device","LightSrcMgr","light-mgr",gtLightMgr);
        currDevice->lightSrcMgr = gpLightMgr[(gpTypeIDptr->second).vectorID];
    }

    //
    // -------- set specific object
    //
    std::string objectID;
    if (gvsParser->getParameter("obj",objectID))
    {
        gpTypeIDptr = getIDptr(gvsParser,"Device","SceneObject","obj",gtSceneObj);
        currDevice->sceneGraph = gpSceneObj[(gpTypeIDptr->second).vectorID];
    }

    //
    // -------- set change values
    //
    std::string objIDname;
    std::string paramName;
    GvsDataType dataType;
    GvsBase* objectPtr;
    void*    paramValue;
    bool     setparamFound = false;

    int num = 0;
    do {
        GvsSetParamList* spl = new GvsSetParamList;
        setparamFound = gvsParser->getParameter("setparam",spl,num);
        if (setparamFound) {
            objIDname  = spl->objectIDname;
            paramName  = spl->paramName;
            dataType   = (spl->param).type;            
            //fprintf(stderr,"Set change value: %s %s %s\n",objIDname.c_str(),paramName.c_str(),GvsDataTypeName[dataType].c_str());

            //std::cerr << objIDname << " " << paramName << std::endl;
            gpTypeIDptr = gpTypeID.find(objIDname);

            if (gpTypeIDptr != gpTypeID.end())
            {
                objectPtr = (gpTypeIDptr->second).gvsObject;
                //objectPtr->Print();

                GvsDataType wantedDataType = objectPtr->GetDataType(paramName);
                dataType = wantedDataType;

                switch (dataType) {
                    case gvsDT_UNKNOWN:
                    case gvsDT_VOID:
                    case gvsDT_STRING:
                        //std::cerr << "parse_device: string\n";
                        break;
                    case gvsDT_INT: {
                        paramValue = (int*)((spl->param).val);
                        break;
                    }
                    case gvsDT_FLOAT: {
                        paramValue = static_cast<float*>((spl->param).val);
                        break;
                    }
                    case gvsDT_DOUBLE: {
                        paramValue = (double*)((spl->param).val);
                        break;
                    }
                    case gvsDT_VEC2: {
                        paramValue = new m4d::vec2( *(static_cast<m4d::vec2*>((spl->param).val)) );
                        break;
                    }
                    case gvsDT_VEC3: {
                        paramValue = new m4d::vec3(*(m4d::vec3*)((spl->param).val));
                        break;
                    }
                    case gvsDT_VEC4: {
                        paramValue = new m4d::vec4(*(m4d::vec4*)((spl->param).val));
                        break;
                    }
                    case gvsDT_MAT2D: {
                        m4d::Matrix<double,2,3>* mptr = (m4d::Matrix<double,2,3>*)((spl->param).val);
                        m4d::Matrix<double,2,3> matrix;
                        for (int i=0; i<2; i++)
                            for (int j=0; j<3; j++)
                                matrix.setCoeff(i,j,mptr->getCoeff(i,j));
                        paramValue = new m4d::Matrix<double,2,3>(matrix);
                        break;
                    }
                    case gvsDT_MAT3D:  {
                        m4d::Matrix<double,3,4>* mptr = (m4d::Matrix<double,3,4>*)((spl->param).val);
                        m4d::Matrix<double,3,4> matrix;
                        for (int i=0; i<3; i++)
                            for (int j=0; j<4; j++)
                                matrix.setCoeff(i,j,mptr->getCoeff(i,j));
                        paramValue = new m4d::Matrix<double,3,4>(matrix);
                        break;
                    }
                    case gvsDT_IVEC2: {
                        paramValue = new m4d::ivec2( *(static_cast<m4d::ivec2*>((spl->param).val)) );
                        //m4d::ivec2* iv = static_cast<m4d::ivec2*>((spl->param).val);
                        //paramValue = new m4d::ivec2( (*iv)[0], (*iv)[1]);
                        break;
                    }
                    case gvsDT_IVEC3: {
                        paramValue = new m4d::ivec3( *(static_cast<m4d::ivec3*>((spl->param).val)) );
                        break;
                    }
                    case gvsDT_IVEC4: {
                        paramValue = new m4d::ivec4( *(static_cast<m4d::ivec4*>((spl->param).val)) );
                        break;
                    }
                }

                currDevice->setChangeObj(objectPtr,objIDname,paramName,dataType, paramValue );
            }
            else
            {
                msg = "ObjectID ";
                msg.append(objIDname);
                msg.append(" not found!\n");
                scheme_error(msg);
            }

        }
        delete spl;
        num++;
    }
    while ( (setparamFound) && (num<gvsParser->getNumParam()) );

    gpDevice.push_back(currDevice);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-device: ID already assigned!");
    }

    GvsTypeID tid = {gtDevice,gpDevice.size()-1,gpDevice[gpDevice.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtDevice"));
    return R;
}
