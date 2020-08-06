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
 *        (csg-object   '(csgtype ,CSGTYPE )
 *                      '(obj1   "obj1" )
 *                      '(obj2   "obj2" )
 *                      '(id    "CSGObjID")
 *        )
 *
 */

#include "Parser/GvsParseScheme.h"
#include "Parser/parse_csgObj.h"
#include "Parser/parse_helper.h"

#include <Obj/SolidObj/GvsSolidObj.h>
#include <Obj/SolidObj/GvsSolidCSGObj.h>
#include <Obj/SolidObj/GvsSolidDifferObj.h>
#include <Obj/SolidObj/GvsSolidIntersecObj.h>
#include <Obj/SolidObj/GvsSolidUnifiedObj.h>
#include <Shader/GvsShader.h>

#include <GvsGlobalDefs.h>

#include "scheme-private.h"

extern std::vector<GvsSceneObj*>    gpSceneObj;
extern std::vector<GvsSolidObj*>    gpSolidObj;
extern std::vector<GvsShader*>      gpShader;

extern std::map<std::string,GvsTypeID>           gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;

//----------------------------------------------------------------------------
//         gvsP_init_solidbox
//----------------------------------------------------------------------------
pointer gvsP_init_csg_obj (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_init_csg_obj..........\n";
#endif

    if (args == sc->NIL) scheme_error("csg-obj: no arguments");
    if (!is_pair(args)) scheme_error("csg-obj: less arguments");

    std::string allowedNames[4] = {"csgType","obj1","obj2","id"};
    GvsParseAllowedNames allowedTypes[4] = {{gp_string_int,1},    // csgType
                                            {gp_string_string,0}, // obj1
                                            {gp_string_string,0}, // obj2
                                            {gp_string_string,0}  // id
                                           };

    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,4);
    args = gvsParser->parse(args);

    GvsSolidObj* child1;
    GvsSolidObj* child2;
    GvsSceneObj* csgObj;
    int          csgType = CSG_Obj;
    std::string CSGTYPE;

    gvsParser->getParameter("csgType",&csgType);

    std::string objID1, objID2;
    gvsParser->getParameter("obj1",objID1);
    gvsParser->getParameter("obj2",objID2);

    //std::cerr << "CSG_Obj1 = " << objID1 << "\tCSG_Obj2 = " << objID2 << std::endl;

    int num = 0;
    if (gpSceneObj.size()>=1) {
        getIDptr(gvsParser,"CSGDifferObject","Object","obj1",gtSceneObj,num);
        child1 = (GvsSolidObj*)gpSceneObj[(gpTypeIDptr->second).vectorID];
    }
    else {
        scheme_error("comp-object: no scene object available!\n");
    }

    if (gpSceneObj.size()>=1) {
        getIDptr(gvsParser,"CSGDifferObject","Object","obj2",gtSceneObj,num);
        child2 = (GvsSolidObj*)gpSceneObj[(gpTypeIDptr->second).vectorID];
    }
    else {
        scheme_error("comp-object: no scene object available!\n");
    }

    if (csgType == CSG_Obj) {
        csgObj  = new GvsSolidDifferObj(child1,child2);
        CSGTYPE = "gpCSGDifferObj";
        std::cerr << "Error: csgtype = " << csgType
                  << ". Invalid CSG_ObjType. CSG_ObjType was set to gpCSGDifferObj."
                  << " Next time choose gpCSGDifferObj = 1, gpCSGUnifiedObj = 2, gpCSGIntersecObj = 3"
                  << std::endl;
        exit(0);
    }

    else if (csgType == CSG_differObj) {
        csgObj = new GvsSolidDifferObj(child1,child2);
        CSGTYPE = "gpCSGDifferObj";
        //      std::cerr << "Error: csgtype = " << CSGTYPE << std::endl;
    }

    else if (csgType == CSG_unifiedObj) {
        csgObj = new GvsSolidUnifiedObj(child1,child2);
        CSGTYPE = "gpCSGUnifiedObj";
        //      std::cerr << "Error: csgtype = " << CSGTYPE << std::endl;
    }

    else if (csgType == CSG_intersecObj) {
        csgObj = new GvsSolidIntersecObj(child1,child2);
        CSGTYPE = "gpCSGIntersecObj";
        //      std::cerr << "Error: csgtype = " << CSGTYPE << std::endl;
    }

    else {
        //      csgObj = new GvsSolidDifferObj(child1,child2);
        std::cerr << "Error: csgtype = " << csgType
                  << ". Invalid CSG_ObjType."
                  << " Next time choose gpCSGDifferObj = 1, gpCSGUnifiedObj = 2, gpCSGIntersecObj = 3"
                  << std::endl;
        CSGTYPE = "gpCSGDifferObj";

        //*/
        exit(0);
    }

    std::string idname;

    gpSceneObj.push_back(csgObj);
    // gpSceneObj[0]->print(cerr);

    idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("csg-obj: ID already assigned!");
    }

    GvsTypeID tid = {gtSceneObj,gpSceneObj.size()-1,gpSceneObj[gpSceneObj.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtSceneObj"));
    return R;
}
