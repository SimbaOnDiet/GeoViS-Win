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
    The projector represents the local observer sitting in the local reference frame
    defined by "locTedObs".
    @verbatim
    (init-projector '(localTetrad "locTedObs")
                   ['(color  #(0.1 0.1 0.1)) ]
                    '(id "proj")
    )@endverbatim
    The default background color is black (0,0,0).


    The local tetrad could also be defined within the initialization of the projector:
    @verbatim
    (init-projector '(pos  #(0.0 10.0 1.5708 0.0 0.0))
                  [ '(e0  #(1.0 0.0 0.0 0.0)) ]
                  [ '(e1  #(0.0 1.0 0.0 0.0)) ]
                  [ '(e2  #(0.0 0.0 1.0 0.0)) ]
                  [ '(e3  #(0.0 0.0 0.0 1.0)) ]
                  [ '(incoords #f) ]
                  ['(color  #(0.1 0.1 0.1)) ]
    )@endverbatim
    If incoords = false, then the local tetrad e0-e3 is given with respect to
    the natural local tetrad.
*/

#include "Parser/parse_projector.h"
#include "Parser/parse_helper.h"

#include "Dev/GvsProjector.h"
#include "Obj/STMotion/GvsLocalTetrad.h"
#include "Obj/STMotion/GvsStMotion.h"

#include "scheme.h"

//#include <Dev/GvsProjector2PI.h>

#include "GvsParseScheme.h"

extern std::vector<GvsLocalTetrad*>    gpLocalTetrad;
extern std::vector<GvsRayGen*>         gpRayGen;
extern std::vector<GvsProjector*>      gpProjector;

extern std::vector<GvsStMotion*>       gpMotion;

extern std::map<std::string,GvsTypeID>           gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;

/**
 * @brief gvsP_init_projector
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_init_projector (scheme *sc, pointer args) {
#ifdef GVS_VERBOSE
    std::cerr << "\n..........gvsP_init_projector..........\n";
#endif
    if (args == sc->NIL) scheme_error("init-projector: no arguments");
    if (!is_pair(args)) scheme_error("init-projector: less arguments");

    GvsProjector*    currProj = NULL;
    GvsLocalTetrad*  locTed = NULL;
    GvsStMotion*     currMotion = NULL;
    //GvsRayGen*       currRayGen = NULL;

    std::string allowedNames[11] = {
        "raygen","localtetrad","color","id","pos","e0","e1","e2","e3","incoords","motion"};
    GvsParseAllowedNames allowedTypes[11] = {{gp_string_string,0}, // raygen
                                             {gp_string_string,0}, // localtetrad
                                             {gp_string_double,3}, // color
                                             {gp_string_string,0}, // id
                                             {gp_string_double,5}, // pos
                                             {gp_string_double,4}, // e0
                                             {gp_string_double,4}, // e1
                                             {gp_string_double,4}, // e2
                                             {gp_string_double,4}, // e3
                                             {gp_string_bool,0},   // incoords
                                             {gp_string_string,0}  // motion
                                            };

    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,11);

    bool haveLocTed = false;
    bool haveMotion = false;


    // Read arguments of projector
    args = gvsParser->parse(args);

    std::string msg;

    // Set ray generator
    std::string rayGenID;
    if (gvsParser->getParameter("raygen",rayGenID)) {
        if (rayGenID=="gtRayGen") {
            currProj = new GvsProjector(gpRayGen[gpRayGen.size()-1]);
        }
        else
        {
            if (gpRayGen.size()>=1)
            {
                gpTypeIDptr = getIDptr(gvsParser,"init-projector","RayGen","raygen",gtRayGen);
                currProj = new GvsProjector(gpRayGen[(gpTypeIDptr->second).vectorID]);
            } else {
                scheme_error("init-projector: no ray generator found!");
            }
        }
    } else {
        if (gpRayGen.empty())  scheme_error("init-projector: RayGen is missing!\n");
        else if (gpRayGen.size()>1) {
            scheme_error("init-projector: RayGen-ID is missing!\n");
        }
        else {
            currProj = new GvsProjector(gpRayGen[0]);
        }
#ifdef GVS_VERBOSE
        printf("Projector: use entered raygen!\n");
#endif
    }

    assert(currProj!=NULL);

    // Determine local tetrad either directly as local tetrad od by position
    double pos[4] = {0.0, 1.0, 0.0, 0.0 };

    std::string locTedID;
    if (gvsParser->getParameter("localtetrad",locTedID)) {
        if (locTedID=="gtLocTed") {
            locTed = gpLocalTetrad[gpLocalTetrad.size()-1];
            haveLocTed = true;
        } else {
            if (gpLocalTetrad.size()>=1) {
                getIDptr(gvsParser,"init-projector","local tetrad","localtetrad",gtLocTed);
                locTed = gpLocalTetrad[(gpTypeIDptr->second).vectorID];
                currProj->setLocalTetrad(locTed);
                haveLocTed = true;
            }
        }
    } else if (gvsParser->getParameter("pos",&pos[0])) {
        currProj->setPosition(m4d::vec4(pos[0],pos[1],pos[2],pos[3]));
        double e0[4] = { 1.0, 0.0, 0.0, 0.0 };
        double e1[4] = { 0.0, 1.0, 0.0, 0.0 };
        double e2[4] = { 0.0, 0.0, 1.0, 0.0 };
        double e3[4] = { 0.0, 0.0, 0.0, 1.0 };
        bool inCoords = false;
        gvsParser->getParameter("e0",&e0[0]);
        gvsParser->getParameter("e1",&e1[0]);
        gvsParser->getParameter("e2",&e2[0]);
        gvsParser->getParameter("e3",&e2[0]);
        gvsParser->getParameter("incoords",inCoords);

        haveLocTed = true;

        std::string locFrame;
        if (gvsParser->getParameter("localFrame",locFrame)) {
            for(int lf=0; lf<3; lf++) {
                //if (GvsLFTypeName[lf]==locFrame) locTed->setInCoords(false,(GvsLFType)lf);
                if (std::string(m4d::stl_nat_tetrad_types[lf])==locFrame)
                    locTed->setInCoords(false,(m4d::enum_nat_tetrad_type)lf);
            }
        }

        currProj->setTetrad(m4d::vec4(e0[0],e0[1],e0[2],e0[3]),
                m4d::vec4(e1[0],e1[1],e1[2],e1[3]),
                m4d::vec4(e2[0],e2[1],e2[2],e2[3]),
                m4d::vec4(e3[0],e3[1],e3[2],e3[3]),
                inCoords);
    }

    // instead of a local tetrad, a motion could also be entered
    std::string motionID;
    if (gvsParser->getParameter("motion",motionID)) {
        gpTypeIDptr = gpTypeID.find(motionID);
        if (gpTypeIDptr != gpTypeID.end()) {
            currMotion = gpMotion[(gpTypeIDptr->second).vectorID];
            haveMotion = true;
        } else {
            msg = "init-projector";
            msg.append(": object with ID \"");
            msg.append(motionID);
            msg.append("\" is no motion!");
            scheme_error(msg);
        }
    }


    if (haveMotion) currProj->setMotion(currMotion);
    // std::cerr << "hier: " << currMotion->getNumPositions() << std::endl;
    // currProj->Print();

    // Set background color
    double color[3] = { 0.0, 0.0, 0.0 };
    if (gvsParser->getParameter("color",&color[0])){
        currProj->setBackgroundColor(GvsColor(color[0],color[1],color[2]));
    }

    if ((!haveLocTed) && (!haveMotion)) {
        msg = "init-projector: ";
        msg.append(locTedID);
        msg.append(" neither a local tetrad nor a motion is given!");
        scheme_error(msg);
    }

    gpProjector.push_back(currProj);

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    } else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-projector: ID already assigned!");
    }

    GvsTypeID tid = {gtProjector,gpProjector.size()-1,gpProjector[gpProjector.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtProjector"));
    return R;
}

