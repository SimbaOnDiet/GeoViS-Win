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
        (init-light-mgr '(ambient #(double double double) )
                        '(light "light-id1" )
                        '(light "light-id2" )
                        '(shadows #t )
                        '(id    "lmgr" )
        )

        (init-light '(type   "PointLight")
                    '(pos    #(double double double double) )
                    '(color  #(double double double) )
                    '(atten  #(double double double) )
                    '(id "pl1")
        )
 */

#include "parse_light.h"
#include "parse_helper.h"

#include "Img/GvsColor.h"
#include "Light/GvsLightSrcMgr.h"
#include "Light/GvsLightSrc.h"
#include "Light/GvsPointLight.h"
#include "Light/GvsFlashLight.h"

#include <GvsGlobalDefs.h>

#include "scheme-private.h"

#include "Parser/GvsParseScheme.h"

extern std::vector<GvsLightSrcMgr*> gpLightMgr;
extern std::vector<GvsLightSrc*>    gpLight;

extern std::map<std::string,GvsTypeID>  gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;



//----------------------------------------------------------------------------
//         gvsP_init_lightmgr
//----------------------------------------------------------------------------
pointer gvsP_init_lightmgr (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n.....gvsP_init_lightmgr.....\n";
#endif
    if (args == sc->NIL) scheme_error("init-light-mgr: no arguments");
    if (!is_pair(args)) scheme_error("init-light-mgr: less arguments");

    std::string allowedNames[5] = {"ambient","diffuse","light","shadows","id"};
    GvsParseAllowedNames allowedTypes[5] = {{gp_string_double,3}, // ambient
                                            {gp_string_double,3}, // diffuse
                                            {gp_string_string,0}, // light
                                            {gp_string_bool,0},   // shadows
                                            {gp_string_string,0}, // id
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,5);
    args = gvsParser->parse(args);

    GvsLightSrcMgr* lightMgr = new GvsLightSrcMgr;

    bool withShadowRays = false;
    gvsParser->getParameter("shadows",withShadowRays);
    lightMgr->setWithShadowRays(withShadowRays);

    double ambient[3] = {1.0, 1.0, 1.0};
    if (gvsParser->getParameter("ambient",ambient)) {
        lightMgr->setAmbientLight(GvsColor(ambient[0],ambient[1],ambient[2]));
    }


    std::string lightID;

    GvsLightSrc* light;
    int num = 0;
    do {
        if (gvsParser->getParameter("light",lightID,num)) {
            if (lightID=="gtLight") {
                light = gpLight[gpLight.size()-1];
                lightMgr->append(light);
            }
            else {
                if (gpLight.size()>=1) {
                    getIDptr(gvsParser,"LightMgr","light","light",gtLight,num);
                    light = gpLight[(gpTypeIDptr->second).vectorID];
                    lightMgr->append(light);
                }
                else
                    scheme_error("LightMgr: no light source available!\n");
            }
        }

        num++;
    } while ( num<gvsParser->getNumParam() );

    gpLightMgr.push_back(lightMgr);

    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtLightMgr"));
    return R;
}

//----------------------------------------------------------------------------
//         gvsP_init_light
//----------------------------------------------------------------------------
pointer gvsP_init_light (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n.....gvsP_init_light.....\n";
#endif
    if (args == sc->NIL) scheme_error("init-light: no arguments");
    if (!is_pair(args)) scheme_error("init-light: less arguments");

    std::string allowedNames[3] = {"type","id"};
    GvsParseAllowedNames allowedTypes[2] = {{gp_string_string,0}, // type
                                            {gp_string_string,0} // id
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,2);
    args = gvsParser->parse(args);

    std::string lightType;
    gvsParser->getParameter("type",lightType);

    if      ( lightType == "PointLight" ) gvsP_init_pointlight(gvsParser);
    else if ( lightType == "FlashLight" ) gvsP_init_flashlight(gvsParser);

    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtLight"));
    return R;

}

//----------------------------------------------------------------------------
//         gvsP_init_pointlight
//----------------------------------------------------------------------------
void gvsP_init_pointlight (GvsParseScheme* gP)
{
    GvsPointLight* pLight = new GvsPointLight;

    double lightPos[4] = {0.0, 6.0, 0.0, 0.0};
    double col[3] = {1.0, 1.0, 1.0};
    //double att[3] = {0.2, 0.0, 0.0};

    gP->setAllowedName("pos",gp_string_double,4);
    gP->setAllowedName("color",gp_string_double,3);
    gP->setAllowedName("atten",gp_string_double,3);

    // Position einlesen
    if (gP->getParameter("pos",lightPos))
        pLight->setPosition(m4d::vec4(lightPos[0],lightPos[1],lightPos[2],lightPos[3]));

    if (gP->getParameter("color",col))
        pLight->setColor(GvsColor(col[0],col[1],col[2]));

    //if (gP->getParameter("atten",att))
    //    pLight->setDistAttenCoeff(m4d::vec3(att[0],att[1],att[2]));


    gpLight.push_back(pLight);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if(gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-light: ID already assigned!");
    }

    GvsTypeID tid = {gtLight,gpLight.size()-1,gpLight[gpLight.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}


//----------------------------------------------------------------------------
//         gvsP_init_flashlight 
//----------------------------------------------------------------------------
void gvsP_init_flashlight (GvsParseScheme* gP)
{
    GvsFlashLight* fLight = new GvsFlashLight;

    double lightPos[5] = {0.0, 6.0, 0.0, 0.0, 0.0};
    double col[3] = {1.0, 1.0, 1.0};
    //double att[3] = {0.2, 0.0, 0.0};

    double lon,loff,synch;
    bool   multi;
    std::string pattern;

    gP->setAllowedName("pos",gp_string_double,5);
    gP->setAllowedName("color",gp_string_double,3);
    gP->setAllowedName("atten",gp_string_double,3);
    gP->setAllowedName("pattern",gp_string_string,0);
    gP->setAllowedName("lengthon",gp_string_double,1);
    gP->setAllowedName("lengthoff",gp_string_double,1);
    gP->setAllowedName("multi",gp_string_bool,0);
    gP->setAllowedName("synchronize",gp_string_double,1);

    // Position einlesen
    if (gP->getParameter("pos",lightPos))
        fLight->setPosition(m4d::vec4(lightPos[0],lightPos[1],lightPos[2],lightPos[3]));

    if (gP->getParameter("color",col))
        fLight->setColor(GvsColor(col[0],col[1],col[2]));

    //if (gP->getParameter("atten",att))
    //    fLight->setDistAttenCoeff(m4d::vec3(att[0],att[1],att[2]));

    if (gP->getParameter("pattern",pattern))
        fLight->setFlashPattern(pattern);

    bool haveOn  = false;
    bool haveOff = false;
    haveOn  = gP->getParameter("lengthon",lon);
    haveOff = gP->getParameter("lengthoff",loff);
    if (haveOn || haveOff)
        fLight->setFlashLengths(lon,loff);

    multi = false;
    if (gP->getParameter("multi",multi))
    {
        if (multi) fLight->setFlashType(gvsFTmultiFlash);
        else fLight->setFlashType(gvsFTsingleFlash);
    }

    if (gP->getParameter("synchronize",synch))
        fLight->setFlashStartTime(synch);

    // fLight->print(cerr);
    gpLight.push_back(fLight);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if(gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-light: ID already assigned!");
    }

    GvsTypeID tid = {gtLight,gpLight.size()-1,gpLight[gpLight.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}
