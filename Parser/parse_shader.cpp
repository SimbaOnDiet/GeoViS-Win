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
        (init-shader '(type  "SurfShader")
                     '(objcolor <texture> )
                     '(ambient   <double> )
                     '(diffuse   <double> )
                     '(id "ID")
        )
 */

#include "parse_shader.h"
#include "parse_helper.h"

#include <GvsGlobalDefs.h>

#include "scheme.h"

#include "Img/GvsColor.h"
#include "Parser/GvsParseScheme.h"
#include "Shader/GvsShader.h"
#include "Shader/Surface/GvsSurfaceShader.h"
#include "Texture/GvsTexture.h"

extern std::vector<GvsShader*>     gpShader;
extern std::vector<GvsTexture*>    gpTexture;

extern std::map<std::string,GvsTypeID>  gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;


//----------------------------------------------------------------------------
//         gvsP_init_shader
//----------------------------------------------------------------------------
pointer gvsP_init_shader (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n.....gvsP_init_shader.....\n";
#endif
    if (args == sc->NIL) scheme_error("init-shader: keine Argumente");
    if (!is_pair(args)) scheme_error("init-shader: zuwenige Argumente");

    std::string allowedNames[3] = {"type","id","color"};
    GvsParseAllowedNames allowedTypes[3] = {{gp_string_string,1}, // type
                                            {gp_string_string,0}, // id
                                            {gp_string_double,3}, // color
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,3);

    // Einlesen der Argumente des Shader
    args = gvsParser->parse(args);

    // Durchsuche Parameter nach dem Typ des Shader
    std::string shaderType;
    gvsParser->getParameter("type",shaderType);

    if (shaderType =="SurfShader") gvsP_init_surfaceShd(gvsParser);

    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtShader"));
    return R;
}



void gvsP_init_surfaceShd (GvsParseScheme* gP)
{
    GvsSurfaceShader* shader = new GvsSurfaceShader;

    std::string apNames[3] = {"objcolor","ambient","diffuse"};
    GvsParseAllowedNames allowedTypes[3] = {{gp_string_string,0}, // objcolor
                                            {gp_string_double,1}, // ambient
                                            {gp_string_double,1}, // diffuse
                                           };
    for (int i=0;i<3;i++) {
        gP->setAllowedName(apNames[i],allowedTypes[i].paramType,allowedTypes[i].paramDim);
    }


    GvsTexture* currTexture = NULL;
    std::string textureID;

    if (gP->getParameter("objcolor",textureID))
    {
        if (textureID=="gtTexture")
            currTexture = gpTexture[gpTexture.size()-1];
        else
        {
            if (gpTexture.size()>=1)
            {
                getIDptr(gP,"init-texture","Texture","texture",gtTexture);
                currTexture = gpTexture[(gpTypeIDptr->second).vectorID];
            }
            else
                scheme_error("init-texture: no texture available!\n");
        }
    } else {
        if (gpTexture.empty()) scheme_error("init-texture: no texture available!\n");
        else if (gpTexture.size()>1)
            scheme_error("init-texture: textur-ID is missing!\n");
        else
            currTexture = gpTexture[(gpTypeIDptr->second).vectorID];
#ifdef GVS_VERBOSE
        printf("init-shader: verwende eingegebene Textur\n");
#endif
    }
    if (currTexture!=NULL) {
        shader->setPrimitiveColor(currTexture);
    }

    // request special parameters
    double ambient, diffuse;
    if (!gP->getParameter("ambient",ambient))  ambient = 1.0;
    if (!gP->getParameter("diffuse",diffuse))  diffuse = 0.0;

    shader->setKambient(ambient);
    shader->setKdiffuse(diffuse);

    gpShader.push_back(shader);
    //shader->Print();

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-shader: ID already assigned!");
    }

    GvsTypeID tid = {gtShader,gpShader.size()-1,gpTexture[gpTexture.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}

