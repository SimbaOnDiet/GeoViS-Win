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
        (solid-box   `(objType ,gpObjTypeInCoords )
                     '(cornerLL  #(-2.0 -2.0 -2.0))
                     '(cornerUR  #( 2.0  2.0  2.0))
                     '(shader    "shader")
                     '(id  "ID")
        )

        (solid-ellipsoid   '(objType     0)
                           '(center   #( 0.0 0.0 0.0))
                           '(axlen    #( 2.0 2.0 2.0))
                           '(shader   "shader")
                           '(id  "ID")
        )


        (solid-background  '(objType     0)
                           '(center   #( 0.0 0.0 0.0))
                           '(axlen    #( 2.0 2.0 2.0))
                           '(shader   "shader")
                           '(id  "ID")
        )


         (solid-cylinder   '(objType     0)
                           '(base     #( 0.0 0.0 0.0))
                           '(top      #( 2.0 2.0 2.0))
                           '(radii    #( 2.0 2.0 2.0 2.0))
                           '(shader   "shader")
                           '(id  "ID")
        )
 */

#include "parse_solidObjects.h"
#include "parse_helper.h"
#include "GvsParseScheme.h"

#include "scheme.h"

#include "GvsGlobalDefs.h"
#include "Obj/SolidObj/GvsSolBox.h"
#include "Obj/SolidObj/GvsSolCylinder.h"
#include "Obj/SolidObj/GvsSolEllipsoid.h"
#include "Obj/SolidObj/GvsSolBackground.h"


#include "Obj/STMotion/GvsLocalTetrad.h"
#include "Obj/STMotion/GvsStMotion.h"

#include <metric/m4dMetric.h>
#include "math/Mat.h"

extern std::vector<GvsLocalTetrad*>    gpLocalTetrad;
extern std::vector<Gvsm4dMetricDummy*> gpMetric;
extern std::vector<GvsShader*>         gpShader;
extern std::vector<GvsSceneObj*>       gpSceneObj;
extern std::vector<GvsStMotion*>       gpMotion;

extern std::map<std::string,GvsTypeID>           gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;


/**
 * @brief gvsP_init_solidbox
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_init_solidbox (scheme *sc, pointer args) {
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_solid-box..........\n";
#endif

    if (args == sc->NIL) scheme_error("solid-box: no arguments");
    if (!is_pair(args)) scheme_error("solid-box: less arguments");

    double cornerLL[3] = { 0.0, 0.0, 0.0 };
    double cornerUR[3] = { 1.0, 1.0, 1.0 };

    GvsSceneObj*  solBox;
    GvsShader*    currShader;
    m4d::Metric*  currMetric;
    int objType = inCoords;

    std::string allowedNames[9] = {"objtype","id","cornerll","cornerur","shader","metric","transform","motion","chart"};
    GvsParseAllowedNames allowedTypes[9] = {{gp_string_int,1},    // objtype
                                            {gp_string_string,0}, // id
                                            {gp_string_double,3}, // cornerLL
                                            {gp_string_double,3}, // cornerUR
                                            {gp_string_string,0}, // shader
                                            {gp_string_string,0}, // metric
                                            {gp_string_matrix,0}, // transform
                                            {gp_string_string,0}, // motion
                                            {gp_string_int,1}     // chart
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,9);

    args = gvsParser->parse(args);
    gvsParser->testParamNames("solid-box");
    gvsParser->getParameter("objtype",&objType);

    gvsParser->getParameter("cornerll",&cornerLL[0]);
    gvsParser->getParameter("cornerur",&cornerUR[0]);

    currShader = readShader("solid-box",gvsParser);
    currMetric = readMetric("solid-box",gvsParser);

    solBox = new GvsSolBox(m4d::vec3(cornerLL[0],cornerLL[1],cornerLL[2]),m4d::vec3(cornerUR[0],cornerUR[1],cornerUR[2]),
            (GvsSurfaceShader*)currShader, currMetric, (GvsObjType)objType);

    int chart = 0;
    if (gvsParser->getParameter("chart",&chart)) {
        solBox->setChart(chart);
    }

    m4d::Matrix<double,3,4> transMat3D;
    if (!gvsParser->getParameter("transform",&transMat3D)) {
        transMat3D.setIdent();
    }
    solBox->transform(transMat3D);


    // read motion: only ConstVelocity is valid
    bool haveMotion = false;
    GvsStMotion* currMotion;
    std::string motionID;

    if (gvsParser->getParameter("motion",motionID)) {        
        gpTypeIDptr = gpTypeID.find(motionID);
        if (gpTypeIDptr != gpTypeID.end()) {
            currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
            haveMotion = true;
        }
        else {
            if ( (gvsParser->getParameter("motion",motionID)) && (gpMotion.size()>0) )  {
                getIDptr(gvsParser,"solid-box","Motion","motion",gtMotion);
                currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
                haveMotion = true;
            }
        }
    }

    if (haveMotion) {
        if (currMotion->getMotionType() == gvsMotionConstVelocity)
            solBox->setMotion(currMotion);
    }

    gpSceneObj.push_back(solBox);
    //gpSceneObj[0]->Print();

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("solid-box: ID is already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}



pointer gvsP_init_solidellipsoid ( scheme *sc, pointer args )
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_solid-ellipsoid..........\n";
#endif

    if (args == sc->NIL) scheme_error("solid-ellipsoid: no arguments");
    if (!is_pair(args)) scheme_error("solid-ellipsoid: less arguments");

    double center[3] = { 0.0, 0.0, 0.0 };
    double axlen[3]  = { 1.0, 1.0, 1.0 };

    GvsSceneObj*   solEllipsoid;
    GvsShader*     currShader;
    m4d::Metric*   currMetric;
    int objType = inCoords;

    std::string allowedNames[9] = {"objtype","id","center","axlen","shader","metric","transform","motion","chart"};
    GvsParseAllowedNames allowedTypes[9] = {{gp_string_int,1},    // objtype
                                            {gp_string_string,0}, // id
                                            {gp_string_double,3}, // center
                                            {gp_string_double,3}, // axlen
                                            {gp_string_string,0}, // shader
                                            {gp_string_string,0}, // metric
                                            {gp_string_matrix,0}, // transform
                                            {gp_string_string,0}, // motion
                                            {gp_string_int,1}     // chart
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,9);
    args = gvsParser->parse(args);

    gvsParser->testParamNames("solid-ellipsoid");
    gvsParser->getParameter("objtype",&objType);
    gvsParser->getParameter("center",&center[0]);
    gvsParser->getParameter("axlen",&axlen[0]);

    currShader = readShader("solid-ellipsoid",gvsParser);
    currMetric = readMetric("solid-ellipsoid",gvsParser);

    solEllipsoid = new GvsSolEllipsoid(m4d::vec3(center[0],center[1],center[2]),
            m4d::vec3(axlen[0],axlen[1],axlen[2]),
            (GvsSurfaceShader*)currShader, currMetric, (GvsObjType)objType);

    int chart = 0;
    if (gvsParser->getParameter("chart",&chart)) {
        solEllipsoid->setChart(chart);
    }
    //solEllipsoid->Print();


    m4d::Matrix<double,3,4> transMat3D;
    if (!gvsParser->getParameter("transform",&transMat3D)) {
        transMat3D.setIdent();
    }
    solEllipsoid->transform(transMat3D);


    // Read motion - only const velocity is permitted
    bool haveMotion = false;
    GvsStMotion* currMotion;
    std::string  motionID;
    if (gvsParser->getParameter("motion",motionID)) {
        if (motionID=="gtMotion") {
            currMotion = gpMotion[gpMotion.size()-1];
            haveMotion = true;
        }
    }
    if ( (gvsParser->getParameter("motion",motionID)) &&
         (gpMotion.size()>0) )
    {
        getIDptr(gvsParser,"local-comp-object","Motion","motion",gtMotion);
        currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
        haveMotion = true;
    }
    if (haveMotion)
        if (currMotion->getMotionType() == gvsMotionConstVelocity)
            solEllipsoid->setMotion(currMotion);

    gpSceneObj.push_back(solEllipsoid);
    // gpSceneObj[0]->print(cerr);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("solid-ellipsoid:  ID already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));

    delete gvsParser;
    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}


//----------------------------------------------------------------------------
//         gvsP_init_solidbackground
//----------------------------------------------------------------------------
pointer gvsP_init_solidbackground ( scheme *sc, pointer args )
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_solid-background..........\n";
#endif

    if (args == sc->NIL) scheme_error("solid-background: no arguments");
    if (!is_pair(args)) scheme_error("solid-background: less arguments");

    GvsSceneObj*     solBackground;
    GvsShader*       currShader;
    m4d::Metric*     currMetric;
    int objType = inCoords;

    std::string allowedNames[5] = {"objtype","id","shader","metric","transform"};
    GvsParseAllowedNames allowedTypes[5] = {{gp_string_int,1},    // objtype
                                            {gp_string_string,0}, // id
                                            {gp_string_string,0}, // shader
                                            {gp_string_string,0}, // metric
                                            {gp_string_matrix,0}  // transform
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,5);

    // Einlesen der Argumente der Metrik
    args = gvsParser->parse(args);
    gvsParser->testParamNames("solid-background");
    gvsParser->getParameter("objtype",&objType);

    currShader = readShader("solid-background",gvsParser);
    currMetric = readMetric("solid-background",gvsParser);

    solBackground = new GvsSolBackground((GvsSurfaceShader*)currShader, currMetric, (GvsObjType)objType);

    // Transformation
    m4d::Matrix<double,3,4> transMat3D;
    bool haveTransformation = gvsParser->getParameter("transform",&transMat3D);
    if (!haveTransformation)
        transMat3D.setIdent();
    solBackground->transform(transMat3D);

    gpSceneObj.push_back(solBackground);
    // gpSceneObj[0]->print(cerr);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("solid-background: ID already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}

/**
 * @brief gvsP_init_solidcylinder
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_init_solidcylinder ( scheme *sc, pointer args )
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_solidcylinder..........\n";
#endif

    if (args == sc->NIL) scheme_error("solid-cylinder: keine Argumente");
    if (!is_pair(args)) scheme_error("solid-clinder: zuwenige Argumente");

    //Hier werden veraenderungen noetig
    double base[3] = { 0.0, 0.0, 0.0 };
    double top[3]  = { 0.0, 0.0, 1.0 };
    double radii[2] = {1.0, 1.0};

    GvsSceneObj*     solCylinder;
    GvsShader*       currShader;
    m4d::Metric*     currMetric;
    int objType = inCoords;

    std::string allowedNames[10] = {"objtype","id","top","base","shader","metric","chart","transform","motion","radii"};
    GvsParseAllowedNames allowedTypes[10] = {{gp_string_int,1},    // objtype
                                             {gp_string_string,0}, // id
                                             {gp_string_double,3}, // top
                                             {gp_string_double,3}, // base
                                             {gp_string_string,0}, // shader
                                             {gp_string_string,0}, // metric
                                             {gp_string_int,1},    // chart
                                             {gp_string_matrix,0}, // transform
                                             {gp_string_string,0}, // motion
                                             {gp_string_double,2}  // radii
                                            };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,10);



    args = gvsParser->parse(args);
    gvsParser->testParamNames("solid-cylinder");
    gvsParser->getParameter("objtype",&objType);
    gvsParser->getParameter("radii",&radii[0]);
    gvsParser->getParameter("top",&top[0]);
    gvsParser->getParameter("base",&base[0]);

    currShader = readShader("solid-cylinder",gvsParser);
    currMetric = readMetric("solid-cylinder",gvsParser);

    solCylinder = new GvsSolCylinder(m4d::vec3(base[0],base[1],base[2]),
            m4d::vec3(top[0],top[1],top[2]),
            m4d::vec2(radii[0],radii[1]),
            (GvsSurfaceShader*)currShader,
            currMetric,
            (GvsObjType)objType);

    int chart = 0;
    if (gvsParser->getParameter("chart",&chart)) {
        solCylinder->setChart(chart);
    }

    // Transformation
    m4d::Matrix<double,3,4> transMat3D;
    if (!gvsParser->getParameter("transform",&transMat3D)) {
        transMat3D.setIdent();
    }
    solCylinder->transform(transMat3D);

    // Read motion: only ConstVelocity is valid.
    GvsStMotion* currMotion = NULL;
    std::string motionID;
    if (gvsParser->getParameter("motion",motionID)) {
        if (motionID=="gtMotion") {
            currMotion = gpMotion[gpMotion.size()-1];
        }
    }
    if ( (gvsParser->getParameter("motion",motionID)) &&
         (gpMotion.size()>0) ) {
        getIDptr(gvsParser,"solid-cylinder","Motion","motion",gtMotion);
        currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
    }
    if (currMotion!=NULL) {
        //currMotion->Print();
        if (currMotion->getMotionType() == gvsMotionConstVelocity)
            solCylinder->setMotion(currMotion);
    }

    gpSceneObj.push_back(solCylinder);
    //solCylinder->Print();

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("solid-cylinder: ID already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}
