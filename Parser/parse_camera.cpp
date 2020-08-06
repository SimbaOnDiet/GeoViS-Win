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

    Parse parameters for a pinhole camera.
    @verbatim
    (init-camera '(type "PinHoleCam")
                 '(dir  #(double double double) )
                 '(vup  #(double double double) )
                 '(fov  #(double double) )
                 '(res  #(int int) )
                 '(filter "string")
                 '(id "string")
    )@endverbatim
    Direction 'dir' and vertical up-vector 'vup' are with respect to the local reference frame
    of the observer - the projector. 'fov' are given in degrees and image resolution 'res' in
    pixels.


    Parse parameters for a 4pi camera. 'angle' defines the main direction and is to be
    given in degrees.
    @verbatim
    (init-camera '(type "4PICam")
                 '(angle  double)
                 '(res  #(int int)
                 '(filter "string")
                 '(id "string")
    )@endverbatim


    Parse parameters for a 2pi camera. The main viewing direction is defined by
    heading and pitch angles (both given in degrees).
    @verbatim
    (init-camera '(type "2PICam")
                 '(heading  double)
                 '(pitch    double)
                 '(res #(int int))
                 '(filter "string")
                 '(id "string")
    )@endverbatim
    Heading and pitch angles are given in degrees.


    Parse parameters for a panorama camera. The field of view (fov) angles have
    to be given in degrees.
    @verbatim
    (init-camera '(type "PanoramaCam")
                 '(dir  #(double double double) )
                 '(vup  #(double double double) )
                 '(fov  #(double double) )
                 '(res  #(int int) )
                 '(filter "string")
                 '(id "string")
    )@endverbatim
    The parameters are the same as for the pinhole camera.


    Camera filters can be "FilterRGB", "FilterRGBpdz", and "FilterRGBjac".
*/


#include "GvsGlobalDefs.h"

#include "Parser/GvsParseScheme.h"
#include "Parser/parse_camera.h"
#include "Parser/parse_helper.h"

#include "Cam/GvsCamera.h"
#include "Cam/GvsPinHoleCam.h"
#include "Cam/GvsPanoramaCam.h"
#include "Cam/Gvs4PICam.h"
#include "Cam/Gvs2PICam.h"

#include "scheme-private.h"


extern std::vector<GvsCamera*>          gpCamera;
extern std::map<std::string,GvsTypeID>  gpTypeID;

/**
 * @brief gvsP_init_camera
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_init_camera (scheme *sc, pointer args) {
#ifdef GVS_VERBOSE
    std::cerr << "\n.....gvsP_init_camera.....\n";
#endif
    if (args == sc->NIL) scheme_error("init-camera: no arguments");
    if (!is_pair(args)) scheme_error("init-camera: less arguments");

    std::string allowedNames[11] = {"type","id","dir","vup","fov","res","filter","param","angle","heading","pitch"};

    GvsParseAllowedNames allowedTypes[11] = {{gp_string_string,0},  // type
                                             {gp_string_string,0},  // id
                                             {gp_string_double,3},  // dir
                                             {gp_string_double,3},  // vup
                                             {gp_string_double,2},  // fov
                                             {gp_string_int,2},     // res
                                             {gp_string_string,0},  // filter
                                             {gp_string_double,1},  // param
                                             {gp_string_double,1},  // angle
                                             {gp_string_double,1},  // heading
                                             {gp_string_double,1}   // pitch
                                           };

    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,11);
    args = gvsParser->parse(args);

    std::string cameraType;
    gvsParser->getParameter("type",cameraType);

    if (cameraType =="PinHoleCam") gvsP_init_pinHoleCam(gvsParser);
    else if (cameraType =="PanoramaCam") gvsP_init_panoramaCam(gvsParser);
    else if (cameraType =="4PICam") gvsP_init_4PICam(gvsParser);
    else if (cameraType =="2PICam") gvsP_init_2PICam(gvsParser);
    else
    {
        std::string msg = cameraType;
        msg.append(": camera is unknown!\n");
        scheme_error(msg);
    }
    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtCamera"));
    return R;
}

/**
 * @brief gvsP_init_pinHoleCam
 * @param gP
 */
void gvsP_init_pinHoleCam (GvsParseScheme* gP) {
    gP->testParamNames("init-camera");

    double dir[3] = {-1.0, 0.0, 0.0};
    double vup[3] = { 0.0, 0.0, 1.0};
    double fov[2] = { 60.0, 48.0 };
    int res[2] = {720,576};
    double param = 0.0;

    gP->getParameter("dir",&dir[0]);
    gP->getParameter("vup",&vup[0]);
    gP->getParameter("fov",&fov[0]);
    gP->getParameter("res",&res[0]);


    // Search filter
    GvsCamFilter filter;
    std::string camFilterName;
    bool haveFilter = gP->getParameter("filter",camFilterName);

    if (haveFilter) {
        int camFilter=-1;
        bool filterFound = false;
        while ((camFilter<=GvsNumCamFilters) && (filterFound==false)) {
            if (camFilterName == GvsCamFilterNames[++camFilter]) filterFound = true;
        }
        if (!filterFound) {
            std::cerr << "This filter does not exist!" << std::endl;
            exit(0);
        }
        filter = GvsCamFilter(camFilter);
    }
    else {
        filter = gvsCamFilterRGB;
    }

    GvsPinHoleCam* phCamera = new GvsPinHoleCam(m4d::vec3(dir[0],dir[1],dir[2]),
            m4d::vec3(vup[0],vup[1],vup[2]),
            m4d::vec2(fov[0],fov[1]),
            m4d::ivec2(res[0],res[1]));
    if (haveFilter) phCamera->setCamFilter(filter);

    if (gP->getParameter("param",&param)) {
        phCamera->setParameter(param);
    }

    gpCamera.push_back(phCamera);
    // phCamera->print(cerr);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-camera: ID already exists!");
    }

    GvsTypeID tid = {gtCamera,gpCamera.size()-1,gpCamera[gpCamera.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}

/**
 * @brief gvsP_init_panoramaCam
 * @param gP
 */
void gvsP_init_panoramaCam( GvsParseScheme* gP ) {
    gP->testParamNames("init-camera");

    double dir[3] = {-1.0, 0.0, 0.0};
    double vup[3] = { 0.0, 0.0, 1.0};
    double fov[2] = { 60.0, 48.0 };
    int res[2] = {720,576};
    double param;

    gP->getParameter("dir",&dir[0]);
    gP->getParameter("vup",&vup[0]);
    gP->getParameter("fov",&fov[0]);
    gP->getParameter("res",&res[0]);


    // Search filter
    GvsCamFilter filter;
    std::string camFilterName;
    bool haveFilter = gP->getParameter("filter",camFilterName);

    if (haveFilter) {
        int camFilter=-1;
        bool filterFound = false;
        while ((camFilter<=GvsNumCamFilters) && (filterFound==false)) {
            if (camFilterName == GvsCamFilterNames[++camFilter]) filterFound = true;
        }
        if (!filterFound) {
            std::cerr << "This filter does not exist!" << std::endl;
            exit(0);
        }
        filter = GvsCamFilter(camFilter);
    }
    else {
        filter = gvsCamFilterRGB;
    }

    GvsPanoramaCam* phCamera = new GvsPanoramaCam(m4d::vec3(dir[0],dir[1],dir[2]),
            m4d::vec3(vup[0],vup[1],vup[2]),
            m4d::vec2(fov[0],fov[1]),
            m4d::ivec2(res[0],res[1]));
    if (haveFilter) phCamera->setCamFilter(filter);
    if (gP->getParameter("param",&param)) {
        phCamera->setParameter(param);
    }

    gpCamera.push_back(phCamera);
    // phCamera->print(cerr);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-camera: ID schon vergeben!");
    }

    GvsTypeID tid = { gtCamera,gpCamera.size()-1,gpCamera[gpCamera.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}

/**
 * @brief gvsP_init_4PICam
 * @param gP
 */
void gvsP_init_4PICam (GvsParseScheme* gP)
{
    gP->testParamNames("init-camera");

    double angle = 0.0;
    int res[2] = {720,576};

    double param = 0.0;
    gP->getParameter("angle",&angle);
    gP->getParameter("res",&res[0]);


    // Search filter
    GvsCamFilter filter;
    std::string camFilterName;
    bool haveFilter = gP->getParameter("filter",camFilterName);

    if (haveFilter) {
        int camFilter=-1;
        bool filterFound = false;
        while ((camFilter<=GvsNumCamFilters) && (filterFound==false)) {
            if (camFilterName == GvsCamFilterNames[++camFilter]) filterFound = true;
        }
        if (!filterFound) {
            std::cerr << "This filter does not exist!" << std::endl;
            exit(0);
        }
        filter = GvsCamFilter(camFilter);
    }
    else {
        filter = gvsCamFilterRGB;
    }

    Gvs4PICam* phCamera = new Gvs4PICam(angle,
                                        m4d::ivec2(res[0],res[1]));
    if (haveFilter) phCamera->setCamFilter(filter);
    if (gP->getParameter("param",&param)) {
        phCamera->setParameter(param);
    }

    gpCamera.push_back(phCamera);
    // phCamera->print(cerr);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if(gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-camera: ID already assigned!");
    }

    GvsTypeID tid = {gtCamera,gpCamera.size()-1,gpCamera[gpCamera.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}

/**
 * @brief gvsP_init_2PICam
 * @param gP
 */
void gvsP_init_2PICam( GvsParseScheme* gP ) {
    gP->testParamNames("init-camera");

    double pitch = 0.0;
    double heading = 0.0;
    int res[2] = {100,100};

    gP->getParameter("heading",&heading);
    gP->getParameter("pitch",&pitch);
    gP->getParameter("res",&res[0]);

    GvsCamFilter filter;
    std::string camFilterName;
    bool haveFilter = gP->getParameter("filter",camFilterName);

    if (haveFilter)
    {
        int camFilter=-1;
        bool filterFound = false;
        do
        {
            if (camFilterName == GvsCamFilterNames[++camFilter]) filterFound = true;
        }
        while ((camFilter<GvsNumCamFilters) && (filterFound==false));
        if (!filterFound) {
            fprintf(stderr,"The filter \"%s\" does not exist!",camFilterName.c_str());
            exit(0);
        }

        filter = GvsCamFilter(camFilter);
    }
    else
    {
        filter = gvsCamFilterRGB;
    }

    Gvs2PICam* phCamera = new Gvs2PICam(heading,pitch,res[0]);
    if (haveFilter) phCamera->setCamFilter(filter);

    gpCamera.push_back(phCamera);
    // phCamera->print(cerr);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-camera: ID already assigned!");
    }

    GvsTypeID tid = {gtCamera,gpCamera.size()-1,gpCamera[gpCamera.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}
