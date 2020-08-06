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

#include "GvsGlobalDefs.h"

#include "Dev/GvsDevice.h"
#include "Obj/GvsSceneObj.h"
#include "Parser/GvsParser.h"

#include "metric/m4dMetric.h"

#include "Utils/GvsLog.h"
extern GvsLog& LOG;


GvsDevice :: GvsDevice() {
    camera      = NULL;
    projector   = NULL;
    metric      = NULL;
    lightSrcMgr = NULL;
    sceneGraph  = NULL;

    listNumChangeObj = 0;
    mChangeObj.clear();
    isManual = false;
}

GvsDevice::~GvsDevice() {
    //fprintf(stderr,"GvsDevice::d'tor\n");
    if (!mChangeObj.empty()) {
        for (unsigned int i=0; i<mChangeObj.size(); i++) {
            if (mChangeObj[i]!=NULL) {
                if (mChangeObj[i]->val != NULL) {
                    if (mChangeObj[i]->type==gvsDT_INT)
                        delete (int*)mChangeObj[i]->val;
                    else if (mChangeObj[i]->type==gvsDT_DOUBLE)
                        delete (double*)mChangeObj[i]->val;
                    else if (mChangeObj[i]->type==gvsDT_VEC3)
                        delete (m4d::vec3*)mChangeObj[i]->val;
                    else if (mChangeObj[i]->type==gvsDT_VEC4)
                        delete (m4d::vec4*)mChangeObj[i]->val;
                    else if (mChangeObj[i]->type==gvsDT_MAT2D)
                        delete (m4d::Matrix<double,2,3>*)mChangeObj[i]->val;
                    else if (mChangeObj[i]->type==gvsDT_MAT3D)
                        delete (m4d::Matrix<double,3,4>*)mChangeObj[i]->val;
                }
                delete mChangeObj[i];
            }
            mChangeObj[i] = NULL;
        }
    }
    mChangeObj.clear();
}


void GvsDevice::setChangeObj( GvsBase* obj, std::string idName, std::string name, GvsDataType tp, void* value )
{
    // std::cerr << "GvsDevice :: setChangeObj" << std::endl;
    GvsChangeObj* chObj = new GvsChangeObj;
    chObj->objectPtr = obj;
    chObj->objectIDname = idName;
    chObj->paramName = name;
    chObj->type = tp;
    chObj->val  = value;

    mChangeObj.push_back(chObj);
    listNumChangeObj = mChangeObj.size();
}


bool GvsDevice::makeChange ( void ) {
    bool changed = true;
    for (int i=0; i<listNumChangeObj; i++) {
        switch (mChangeObj[i]->type) {
            case gvsDT_INT: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(int*)mChangeObj[i]->val);
                break;
            }
            case gvsDT_DOUBLE: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(double*)mChangeObj[i]->val);
                break;
            }
            case gvsDT_VEC2: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(static_cast<m4d::vec2*>(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_VEC3: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(static_cast<m4d::vec3*>(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_VEC4: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(static_cast<m4d::vec4*>(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_MAT2D: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(m4d::Matrix<double,2,3>*)mChangeObj[i]->val);
                break;
            }
            case gvsDT_MAT3D: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(m4d::Matrix<double,3,4>*)mChangeObj[i]->val);
                break;
            }
            case gvsDT_IVEC2: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(static_cast<m4d::ivec2*>(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_IVEC3: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(static_cast<m4d::ivec3*>(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_IVEC4: {
                (mChangeObj[i]->objectPtr)->SetParam(mChangeObj[i]->paramName,*(static_cast<m4d::ivec4*>(mChangeObj[i]->val)));
                break;
            }
            default: {
                changed = false;
                break;
            }
        } // end switch
    }


    // THE FOLLOWING IS NOT ALLOWED WHEN SETTING TETRAD VECTORS MANUALLY !!!!
    // In that case, you have to call
    //   (setparam ("tedID" "calc" 0))
    // to manually adjust tetrad when a metric parameter has changed.

    // If the projector is described by a motion, then it is not allowed
    // to adjust the local tetrad !!
    if (!isManual) {
        GvsStMotion* motion = projector->getMotion();
        GvsLocalTetrad* lt = projector->getLocalTetrad();
        if (motion==NULL && lt!=NULL) {
            lt->adjustTetrad();
        }
    }
    return changed;
}


void GvsDevice :: clear ( void ) {
    camera      = NULL;
    projector   = NULL;
    metric      = NULL;
    lightSrcMgr = NULL;
    sceneGraph  = NULL;
}


void GvsDevice :: setManual ( bool manual ) {
    isManual = manual;
}

void GvsDevice::Print( FILE* fptr ) {
    LOG.printf("Device consists of the following components:\n");
    LOG.printf("--------------------------------------------\n");
    if (camera!=NULL) {
        camera->Print(fptr);
    }
    if (projector!=NULL) {
        projector->Print(fptr);
    }
    if (metric!=NULL) {
        metric->print(fptr);
    }
    if (lightSrcMgr!=NULL) {
        lightSrcMgr->Print(fptr);
    }
    if (sceneGraph!=NULL) {
        sceneGraph->Print(fptr);
    }

    LOG.printf("\nChange:\n");
    for(unsigned int i=0; i<mChangeObj.size(); i++) {
        LOG.printf("\tIDname    [%d]:  %s\n",i,mChangeObj[i]->objectIDname.c_str());
        LOG.printf("\tParamName [%d]:  %s\n",i,mChangeObj[i]->paramName.c_str());
        LOG.printf("\tParamType [%d]:  %s\n",i, GvsDataTypeName[mChangeObj[i]->type].c_str());
        LOG.printf("\tParamValue[%d]:  ",i);
        switch (mChangeObj[i]->type) {
            case gvsDT_VOID:
            case gvsDT_STRING:
            case gvsDT_UNKNOWN: {
                break;
            }
            case gvsDT_INT: {
                LOG.printf("%d\n",(int)(*(int*)(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_FLOAT: {
                LOG.printf("%f\n",(float)(*(float*)(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_DOUBLE: {
                LOG.printf("%f\n",(double)(*(double*)(mChangeObj[i]->val)));
                break;
            }
            case gvsDT_VEC2: {
                m4d::vec2* v = (m4d::vec2*)(mChangeObj[i]->val);
                // TODO: move output to logfile
                v->printS(fptr);
                break;
            }
            case gvsDT_VEC3: {
                m4d::vec3* v = (m4d::vec3*)(mChangeObj[i]->val);
                v->printS(fptr);
                break;
            }
            case gvsDT_VEC4: {
                m4d::vec4* v = (m4d::vec4*)(mChangeObj[i]->val);
                v->printS(fptr);
                break;
            }
            case gvsDT_MAT2D: {
                m4d::Matrix<double,2,3>* mat = (m4d::Matrix<double,2,3>*)(mChangeObj[i]->val);
                mat->printS(fptr);
                break;
            }
            case gvsDT_MAT3D: {
                m4d::Matrix<double,3,4>* mat = (m4d::Matrix<double,3,4>*)(mChangeObj[i]->val);
                mat->printS(fptr);
                break;
            }
            case gvsDT_IVEC2: {
                m4d::ivec2* v = static_cast<m4d::ivec2*>(mChangeObj[i]->val);
                v->printS(fptr);
                break;
            }
            case gvsDT_IVEC3: {
                m4d::ivec3* v = static_cast<m4d::ivec3*>(mChangeObj[i]->val);
                v->printS(fptr);
                break;
            }
            case gvsDT_IVEC4: {
                m4d::ivec4* v = static_cast<m4d::ivec4*>(mChangeObj[i]->val);
                v->printS(fptr);
                break;
            }
        }
    }
    LOG.printf("End of Change.\n");
}
