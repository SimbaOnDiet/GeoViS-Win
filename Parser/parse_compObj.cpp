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
    @brief  Parse compound object- source
    @file   parse_compObj.cpp

    @verbatim
    (comp-obj '(obj "objID1")
              '(obj "objID2")
              ...
              `(transform  ...)
              '(id "string")
    )@endverbatim
    Note that the transformation only acts on already included objects. Those objects
    that are added later by 'add-object' will not be transformed!


    @verbatim
    (local-comp-obj '(obj "objID1")
                    '(obj "objID2")
                    '(localtetrad "locTedID")
                    '(motion "motionID")
                    '(id "string")
    )@endverbatim
    A local compount object is being described either by a motion OR by a
    local tetrad. The latter one can be represented in coordinates as well
    as with respect to a natural local tetrad.


    @verbatim
    (add-object  '(add-to "compObjID")
                 '(obj  "objID")
    )@endverbatim
*/

#include <Parser/parse_compObj.h>
#include <Parser/parse_locTed.h>
#include <Parser/parse_helper.h>
#include "GvsParseScheme.h"

#include <Obj/STMotion/GvsStMotion.h>
#include <Obj/Comp/GvsLocalCompObj.h>
#include <Obj/Comp/GvsCompoundObj.h>

#include "scheme-private.h"


extern std::vector<GvsLocalTetrad*>     gpLocalTetrad;
extern std::vector<GvsLocalCompObj*>    gpLocalCompObj;
extern std::vector<GvsSceneObj*>        gpSceneObj;

extern std::vector<GvsStMotion*>        gpMotion;

extern std::map<std::string,GvsTypeID>           gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;


/**
 * @brief gvsP_compound_obj
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_compound_obj (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_compound_obj..........\n";
#endif
    if (args == sc->NIL) scheme_error("comp-object: no arguments");
    if (!is_pair(args)) scheme_error("comp-object: less arguments");

    std::string allowedNames[3] = {"id","obj","transform"};
    GvsParseAllowedNames allowedTypes[3] = {{gp_string_string,0},   // id
                                            {gp_string_string,0},  // object
                                            {gp_string_matrix,0}   // transform
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,3);
    args = gvsParser->parse(args);
    gvsParser->testParamNames("comp-object");

    GvsCompoundObj* compObj = new GvsCompoundObj;

    std::string objID;
    GvsSceneObj* object;

    // if (gpSceneObj.empty()) { scheme_error("CompObj: Szene-Objekt fehlt\n"); }
    int num = 0;
    do {
        if (gvsParser->getParameter("obj",objID,num)) {
            if (objID=="gtSceneObj") {
                object = gpSceneObj[gpSceneObj.size()-1];
                compObj->Add(object);
            }
            else {
                if (gpSceneObj.size()>=1) {
                    getIDptr(gvsParser,"CompObj","Object","obj",gtSceneObj,num);
                    object = gpSceneObj[(gpTypeIDptr->second).vectorID];
                    compObj->Add(object);
                }
                else {
                    scheme_error("comp-object: no scene object available!\n");
                }
            }
        }
        num++;
    }
    while ( num<gvsParser->getNumParam() );

    if (compObj->getNumObjs() == 0) {
        std::cerr << "Warning: no scene object read!" << std::endl;
    }

    m4d::Matrix<double,3,4> transMat3D;
    if (!gvsParser->getParameter("transform",&transMat3D)) {
        transMat3D.setIdent();
    }
    compObj->transform(transMat3D);

    gpSceneObj.push_back(compObj);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("comp-object: ID already assigned!");
    }

    GvsTypeID tid = { gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}


/**
 * @brief gvsP_local_comp_obj
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_local_comp_obj (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_local_comp_obj..........\n";
#endif
    if (args == sc->NIL) scheme_error("local-comp-object: no arguments");
    if (!is_pair(args)) scheme_error("local-comp-object: less arguments");


    std::string allowedNames[4] = {"id","localtetrad","obj","motion"};
    GvsParseAllowedNames allowedTypes[4] = {{gp_string_string,0}, // id
                                            {gp_string_string,0}, // localtetrad
                                            {gp_string_string,0}, // obj
                                            {gp_string_string,0}  // motion
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,4);
    args = gvsParser->parse(args);
    gvsParser->testParamNames("local-comp-obj");

    bool haveLocalTetrad = false;
    bool haveMotion      = false;

    GvsLocalCompObj* locCompObj = new GvsLocalCompObj;
    GvsStMotion*     currMotion = NULL;
    GvsLocalTetrad*  currLocTed = NULL;


    // A local compound object is described either by a local tetrad (static object)
    // or by a motion. Both possibilities have to be checked.

    // -- search for local tetrad
    std::string locTedID;
    if (gvsParser->getParameter("localtetrad",locTedID)) {
        if (locTedID=="gtLocTed") {
            currLocTed = gpLocalTetrad[gpLocalTetrad.size()-1];
            haveLocalTetrad = true;
        }
        else {
            if (gpLocalTetrad.size()>0) {
                gpTypeIDptr = getIDptr(gvsParser,"local-comp-object","LocalTetrad","localtetrad",gtLocTed);
                currLocTed = gpLocalTetrad[(gpTypeIDptr->second).vectorID];
                haveLocalTetrad = true;
            }
        }
    }


    // -- search for motion
    std::string motionID;
    if (gvsParser->getParameter("motion",motionID)) {
        if (motionID=="gtMotion") {
            currMotion = gpMotion[gpMotion.size()-1];
            haveMotion = true;
        }
        else {
            if (gpMotion.size()>=1) {
                gpTypeIDptr = getIDptr(gvsParser,"local-comp-obj","Motion","motion",gtMotion);
                currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
                haveMotion = true;
            }
        }
    }


    // Teste, ob Bewegung und oder Lokale Tetrade eingegeben wurde
    if (!haveLocalTetrad && !haveMotion)
        scheme_error("local-comp-object: Lokale Tetrade oder Bewegung eingeben!");

    if (haveLocalTetrad && haveMotion)
        scheme_error("local-comp-object: Entweder lokale Tetrade oder Bewegung eingeben!");


    if (haveLocalTetrad)
    {
#ifdef GVS_VERBOSE
        std::cerr << "\n   lokale Tetrade gefunden...\n";
#endif
        currLocTed->adjustTetrad(); // ------------ !! noch anpassen mit Abfragen --------------------
        locCompObj->setLocalTetrad(currLocTed);
    }
    else if (haveMotion)
    {
#ifdef GVS_VERBOSE
        std::cerr << "\n   Bewegung gefunden...\n";
#endif
        locCompObj->setMotion(currMotion);
        currMotion->Print();
    }

    // read all objects
    std::string objID;
    bool objFound = false;
    GvsSceneObj* object;

    int num = 0;
    do {
        objFound = false;

        //   std::cerr << "objID: " << objID << std::endl;
        if (gvsParser->getParameter("obj",objID,num))
        {
            if (objID=="gtSceneObj")
            {
                object = gpSceneObj[gpSceneObj.size()-1];
                objFound = true;
            }
        }
        if ( (gvsParser->getParameter("obj",objID,num)) &&
             (gpSceneObj.size()>0) )
        {
            gpTypeIDptr = getIDptr(gvsParser,"local-comp-object","Object","obj",gtSceneObj,num);
            object = gpSceneObj[(gpTypeIDptr->second).vectorID];
            //object->Print();
            objFound = true;
        }

        if (objFound) {
            locCompObj->Add(object);
        }

        num++;
    }
    while ( (objFound) && (num<gvsParser->getNumParam()));

    if (locCompObj->getNumObjs()==0)
        scheme_error("local-comp-object: no object available!");
    locCompObj->calcSTBoundBoxComplete();

    //locCompObj->Print();

    gpSceneObj.push_back(locCompObj);

    // In der Objekt-Map einen Eintrag schreiben, wo sich das LocalCompObj im gpSceneObj-Vektor befindet
    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("local-comp-obj: ID already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;
    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}


/**
 * @brief gvsP_add_object
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_add_object ( scheme *sc, pointer args )
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_add_object..........\n";
#endif
    if (args == sc->NIL) scheme_error("add-object: no arguments");
    if (!is_pair(args)) scheme_error("add-object: less arguments");

    std::string allowedNames[2] = {"add-to","obj"};
    GvsParseAllowedNames allowedTypes[2] = {{gp_string_string,0}, // add-to
                                            {gp_string_string,0}  // object
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,2);
    args = gvsParser->parse(args);
    gvsParser->testParamNames("add-object");

    std::string msg;
    std::string idName;
    if (!gvsParser->getParameter("add-to",idName)) {
        scheme_error("add-object: add-to ID is missing!");
    }

    gpTypeIDptr = gpTypeID.find(idName);
    if (gpTypeIDptr == gpTypeID.end()) {
        msg = "ID: ";
        msg.append(idName);
        msg.append(" is not available!");
        scheme_error(msg);
    }

    // ... wenn ja
    //int numParams = gvsParser->getNumParam();
    //std::string paramName;
    //int dim;

    GvsBase* gvsObject = (gpTypeIDptr->second).gvsObject;

    if ((gpTypeIDptr->second).gvsType != gtSceneObj)
    {
        msg = "add-object: ";
        msg.append(idName);
        msg.append(" is no scene object!");
        scheme_error(msg);
    }

    std::string objID;
    bool objFound = false;
    GvsSceneObj* object = NULL;

    int num = 0;
    do {
        objFound = false;

        //   std::cerr << "objID: " << objID << std::endl;
        if (gvsParser->getParameter("obj",objID,num)) {
            if (objID=="gtSceneObj") {
                object = gpSceneObj[gpSceneObj.size()-1];
                objFound = true;
            }
        }
        if ( (gvsParser->getParameter("obj",objID,num)) &&
             (gpSceneObj.size()>0) )
        {
            getIDptr(gvsParser,"add-object","Object","obj",gtSceneObj,num);
            object = gpSceneObj[(gpTypeIDptr->second).vectorID];
            //	object->print(cerr);
            objFound = true;
        }

        if (objFound) {
            // std::cerr << "addiere Objekt\n";
            gvsObject->Add(object);
        }

        num++;
    }
    while ( (objFound) && (num<gvsParser->getNumParam()));

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}
