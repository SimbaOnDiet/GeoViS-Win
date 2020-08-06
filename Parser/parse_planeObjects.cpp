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
        (plane-triangle  `(objtype ,gpObjTypeInCoords)
                         '(p1  #(double double double))
                         '(p2  #(double double double))
                         '(p3  #(double double double))
                         '(shader "shaderID")
                       [ '(metric "metricID") ]
                       [ '(motion "motionID") ]
                       [ '(chart  integer) ]
                         '(id    "triangleID" )
        )


        (plane-ring `(objtype ,gpObjTypeInCoords)
                    '(center  #(double double double))
                    '(normal  #(double double double))
                    '(rout  double)
                    '(rin   double)
                    '(shader "shaderID")
                  [ '(metric "metricID") ]
                  [ '(motion "motionID") ]
                  [ '(chart  integer) ]
                    '(id    "ringID" )
        )
  */


#include "parse_planeObjects.h"
#include "parse_helper.h"
#include "GvsParseScheme.h"

#include "GvsGlobalDefs.h"
#include "Obj/PlanarObj/GvsRing.h"
#include "Obj/PlanarObj/GvsTriangle.h"

#include "Obj/STMotion/GvsLocalTetrad.h"
#include "Obj/STMotion/GvsStMotion.h"

#include "metric/m4dMetric.h"
#include "math/Mat.h"

#include "scheme-private.h"

extern std::vector<GvsLocalTetrad*>    gpLocalTetrad;
extern std::vector<Gvsm4dMetricDummy*> gpMetric;
extern std::vector<GvsShader*>         gpShader;
extern std::vector<GvsSceneObj*>       gpSceneObj;
extern std::vector<GvsStMotion*>       gpMotion;

extern std::map<std::string,GvsTypeID>           gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;



pointer gvsP_init_planetriangle (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_plane-triangle..........\n";
#endif

    if (args == sc->NIL) scheme_error("plane-triangle: no arguments");
    if (!is_pair(args)) scheme_error("plane-triangle: less arguments");


    GvsSceneObj*     planeTriangle;
    GvsShader*       currShader;
    m4d::Metric*     currMetric;
    int objType = inCoords;

    std::string allowedNames[9] = {"objtype","id","p1","p2","p3","shader","metric","motion","chart"};
    GvsParseAllowedNames allowedTypes[9] = {{gp_string_int,1},    // objtype
                                             {gp_string_string,0}, // id
                                             {gp_string_double,3}, // p1
                                             {gp_string_double,3}, // p2
                                             {gp_string_double,3}, // p3
                                             {gp_string_string,0}, // shader
                                             {gp_string_string,0}, // metric
                                             {gp_string_string,0}, // motion
                                             {gp_string_int,1}     // chart
                                            };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,9);

    args = gvsParser->parse(args);
    gvsParser->testParamNames("plane-triangle");
    gvsParser->getParameter("objtype",&objType);

    double p1[3],p2[3],p3[3];
    gvsParser->getParameter("p1",&p1[0]);
    gvsParser->getParameter("p2",&p2[0]);
    gvsParser->getParameter("p3",&p3[0]);

    currShader = readShader("plane-triangle",gvsParser);
    currMetric = readMetric("plane-triangle",gvsParser);

    planeTriangle = new GvsTriangle(m4d::vec3(p1[0],p1[1],p1[2]),
            m4d::vec3(p2[0],p2[1],p2[2]),
            m4d::vec3(p3[0],p3[1],p3[2]),
            (GvsSurfaceShader*)currShader, currMetric, (GvsObjType)objType);

    int chart = 0;
    if (gvsParser->getParameter("chart",&chart)) {
        planeTriangle->setChart(chart);
    }

    bool haveMotion = false;
    GvsStMotion* currMotion;
    std::string motionID,msg;
    if (gvsParser->getParameter("motion",motionID)) {

        gpTypeIDptr = gpTypeID.find(motionID);
        if (gpTypeIDptr != gpTypeID.end()) {
            currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
            haveMotion = true;
        } else {
            msg = "plane-triangle";
            msg.append(": object with ID \"");
            msg.append(motionID);
            msg.append("\" is no motion!");
            scheme_error(msg);
        }
    }

    if (haveMotion)  {
        if (currMotion->getMotionType() == gvsMotionConstVelocity)
            planeTriangle->setMotion(currMotion);
    }

    gpSceneObj.push_back(planeTriangle);
    // gpSceneObj[0]->print(cerr);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("plane-triangle: ID is already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}



pointer gvsP_init_planering (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_plane-ring..........\n";
#endif

    if (args == sc->NIL) scheme_error("plane-ring: no arguments");
    if (!is_pair(args)) scheme_error("plane-ring: less arguments");


    GvsSceneObj*     planeRing;
    GvsShader*       currShader;
    m4d::Metric*     currMetric;
    int objType = inCoords;

    std::string allowedNames[10] = {"objtype","id","center","normal","rout","rin","shader","metric","motion","chart"};
    GvsParseAllowedNames allowedTypes[10] = {{gp_string_int,1},    // objtype
                                             {gp_string_string,0}, // id
                                             {gp_string_double,3}, // center
                                             {gp_string_double,3}, // normal
                                             {gp_string_double,1}, // rout
                                             {gp_string_double,1}, // rin
                                             {gp_string_string,0}, // shader
                                             {gp_string_string,0}, // metric
                                             {gp_string_string,0}, // motion
                                             {gp_string_int,1}     // chart
                                            };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,10);

    args = gvsParser->parse(args);
    gvsParser->testParamNames("plane-ring");
    gvsParser->getParameter("objtype",&objType);

    double center[3],normal[3],rout,rin;
    gvsParser->getParameter("center",&center[0]);
    gvsParser->getParameter("normal",&normal[0]);
    gvsParser->getParameter("rout",&rout);
    gvsParser->getParameter("rin",&rin);

    currShader = readShader("plane-ring",gvsParser);
    currMetric = readMetric("plane-ring",gvsParser);

    planeRing = new GvsRing(m4d::vec3(center[0],center[1],center[2]),
            m4d::vec3(normal[0],normal[1],normal[2]),
            rout,rin,
            (GvsSurfaceShader*)currShader, currMetric, (GvsObjType)objType);

    int chart = 0;
    if (gvsParser->getParameter("chart",&chart)) {
        planeRing->setChart(chart);
    }


    bool haveMotion = false;
    GvsStMotion* currMotion;
    std::string motionID,msg;
    if (gvsParser->getParameter("motion",motionID)) {

        gpTypeIDptr = gpTypeID.find(motionID);
        if (gpTypeIDptr != gpTypeID.end()) {
            currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
            haveMotion = true;
        } else {
            msg = "plane-ring";
            msg.append(": object with ID \"");
            msg.append(motionID);
            msg.append("\" is no motion!");
            scheme_error(msg);
        }
    }

    if (haveMotion) {
        if (currMotion->getMotionType() == gvsMotionConstVelocity)
            planeRing->setMotion(currMotion);
    }

    gpSceneObj.push_back(planeRing);
    // gpSceneObj[0]->print(cerr);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("plane-ring: ID is already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}


