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
    Initializing a solver can be done in the following way. Arguments in brackets [ ] are optional.
    @verbatim
    (init-solver '(type "SolverName")
               [ '(geodType "lightlike") ]
               [ '(geodDir "forward")    ]
               [ '(step_ctrl #t)         ]
               [ '(step_size <double>)   ]
               [ '(max_step <double>)    ]
               [ '(eps_abs <double>)     ]
               [ '(eps_rel <double>)     ]
               [ '(id "solver")          ]
    )@endverbatim

    - The solver type depends on the Motion4D library. Possible values are 'GSL_RK_Cash-Karp',...
    - The geodesic type (geodType) can be either 'lightlike' or 'timelike'.
    - The direction can only be 'forward' or 'backward'. If the solver is used for raytracing, then
      the direction is automatically set to 'backward'.
*/

#include "Parser/parse_solver.h"
#include "Parser/parse_helper.h"
#include "GvsParseScheme.h"

#include <metric/m4dMetric.h>
#include <motion/m4dMotionDatabase.h>
#include "scheme.h"

#include "Utils/GvsGeodSolver.h"

extern std::vector<GvsGeodSolver*>  gpSolver;
extern std::vector<Gvsm4dMetricDummy*>  gpMetric;

extern std::map<std::string,GvsTypeID>           gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;

/**
 * @brief gvsP_init_solver
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_init_solver (scheme *sc, pointer args) {
#ifdef GVS_VERBOSE
    std::cerr << "\n.....gvsP_init_solver.....\n";
#endif

    if (args == sc->NIL) scheme_error("init-solver: no arguments");
    if (!is_pair(args)) scheme_error("init-solver: less arguments");

    std::string allowedNames[12] = {
        "type","metric","geodtype","geoddir","step_ctrl","step_size","max_step",
        "eps_abs","eps_rel","boundboxll","boundboxur","id"};

    GvsParseAllowedNames allowedTypes[12] = {{gp_string_string,0},  // type
                                             {gp_string_string,0},  // metric
                                             {gp_string_string,0},  // geodtype
                                             {gp_string_string,0},  // geoddir
                                             {gp_string_bool,0},    // step_ctrl
                                             {gp_string_double,1},  // step_size
                                             {gp_string_double,1},  // max_step
                                             {gp_string_double,1},  // eps_abs
                                             {gp_string_double,1},  // eps_rel
                                             {gp_string_double,4},  // boundBoxLL
                                             {gp_string_double,4},  // boundBoxUR
                                             {gp_string_string,0}   // id
                                            };

    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,12);
    args = gvsParser->parse(args);
    gvsParser->testParamNames("init-solver");

    GvsGeodSolver* currSolver = NULL;
    m4d::Metric*   currMetric = NULL;

    currMetric = readMetric("init-solver",gvsParser);
    if (currMetric==NULL) {
        exit(-1);
    }

    std::string solverName;
    if (!gvsParser->getParameter("type",solverName)) {
        scheme_error("init-solver: type is missing!");
    }

    m4d::IntegratorDatabase* IntDB = new m4d::IntegratorDatabase;
    m4d::enum_integrator solverID = IntDB->getIntegratorNr(solverName);
    if (solverID==m4d::gsUnknown) scheme_error("Solver is not in the m4d database!");
    delete IntDB;

    // Initialize solver
#ifdef GVS_VERBOSE
    printf("\n-->Initialize solver...\n");
#endif

    if (solverID!=m4d::gsUnknown){
        currSolver = new GvsGeodSolver(currMetric, solverID);
    }

    // Set geodesic type
    std::string geodType;
    if (gvsParser->getParameter("geodtype",geodType)) {
        if (geodType == "lightlike") currSolver->setGeodType(m4d::enum_geodesic_lightlike);
        else if (geodType == "timelike") currSolver->setGeodType(m4d::enum_geodesic_timelike);
        else if (geodType == "spacelike") currSolver->setGeodType(m4d::enum_geodesic_spacelike);  // not tested
    } else {
        currSolver->setGeodType(m4d::enum_geodesic_lightlike);
    }

    // Set geodesic direction
    std::string geodDir;
    if (gvsParser->getParameter("geoddir",geodDir)) {
        if (geodDir == "forward") currSolver->setTimeDir(m4d::enum_time_forward);
        else if (geodDir == "backward") currSolver->setTimeDir(m4d::enum_time_backward);
        else scheme_error("init-solver: wrong time direction!");
    }

    bool step_control;
    if (!gvsParser->getParameter("step_ctrl",step_control)) step_control = true;
    currSolver->setStepSizeControl(step_control);

    double step_size;
    if (!gvsParser->getParameter("step_size",step_size)) step_size = 0.01;
    currSolver->setStepsize(step_size);

    double maxStep;
    if (!gvsParser->getParameter("max_step",maxStep)) maxStep = DEF_MAX_STEPSIZE;
    currSolver->setMaxStepsize(maxStep);


    double eps_abs,eps_rel;
    if (!gvsParser->getParameter("eps_abs",eps_abs)) eps_abs = 1.0e-4;
    if (!gvsParser->getParameter("eps_rel",eps_rel)) eps_rel = 0.0;
    currSolver->setEpsilons(eps_abs,eps_rel);

    // Set bounding box
    double boundBoxLL[4] = {-DBL_MAX,-50.0,-50.0,-50.0};
    double boundBoxUR[4] = { DBL_MAX, 50.0, 50.0, 50.0};

    bool haveBoxLL = gvsParser->getParameter("boundboxll",&boundBoxLL[0]);
    bool haveBoxUR = gvsParser->getParameter("boundboxur",&boundBoxUR[0]);
    if (haveBoxLL || haveBoxUR) {
        currSolver->setBoundingBox(boundBoxLL,boundBoxUR);
    }

    gpSolver.push_back(currSolver);

#ifdef GVS_VERBOSE
    currSolver->Print();
    printf("\n");
#endif

    std::string idname = "unknown";
    if (!gvsParser->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-solver: ID already assigned!");
    }

    GvsTypeID tid = {gtGeodSolver,gpSolver.size()-1,gpSolver[gpSolver.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));    
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtGeodSolver"));
    return R;
}
