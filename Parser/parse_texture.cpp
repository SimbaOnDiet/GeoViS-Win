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
        (init-texture '(type  "UniTex")
                      '(color #(double double double) )
                    [ '(colorgray 0.3) ]
        )

        (init-texture '(type  "CheckerT2D")
                      '(texture "texID")
                      '(texture "texID")
                      `(transform ,(scale-obj #(2.0 2.0)))
                      '(id "tex1")
        )

        (init-texture '(type  "ChequeredT2D")
                      '(texture "texID")
                      '(texture "texID")
                      '(width  0.1)
                      `(transform ,(scale-obj #(2.0 2.0)))
                      '(id "tex1")
        )

        (init-texture '(type "Image")
                      '(file "filename")
                      `(transform ,(scale-obj #(2.0 2.0)))
                      '(id "imageTex")
        )
 */


#include "parse_texture.h"
#include "parse_helper.h"

#include <GvsGlobalDefs.h>

#include "Img/GvsColor.h"
#include "Img/GvsChannelImg2D.h"
#include "Img/GvsPicIOEnvelope.h"
#include "Parser/GvsParseScheme.h"
#include "Texture/GvsTexture.h"
#include "Texture/GvsUniTex.h"
#include <Texture/GvsImg2DSampler.h>

#include "scheme.h"

#include <Texture/GvsCheckerT2D.h>
//#include <Texture/2D/GvsChequeredT2D.h>
//#include <Texture/2D/GvsChequeredTimeT2D.h>


#include "math/Mat.h"


extern std::vector<GvsTexture*>        gpTexture;
extern std::vector<GvsChannelImg2D*>   gpChannelImg2D;

extern std::map<std::string,GvsTypeID>  gpTypeID;
extern std::map<std::string,GvsTypeID>::iterator gpTypeIDptr;

/**
 * @brief gvsP_init_texture
 * @param sc
 * @param args
 * @return
 */
pointer gvsP_init_texture (scheme *sc, pointer args)
{
#ifdef GVS_VERBOSE
    std::cerr << "\n.....gvsP_init_texture.....\n";
#endif
    if (args == sc->NIL) scheme_error("init-texture: no arguments");
    if (!is_pair(args)) scheme_error("init-texture: less arguments");

    std::string allowedNames[5] = {"type","id","color","colorgray","transform"};
    GvsParseAllowedNames allowedTypes[5] = {{gp_string_string,0}, // type
                                            {gp_string_string,0}, // id
                                            {gp_string_double,3}, // color
                                            {gp_string_double,1}, // colorgray
                                            {gp_string_matrix,0}, // transform
                                           };
    GvsParseScheme* gvsParser = new GvsParseScheme(sc,allowedNames,allowedTypes,5);

    // Einlesen der Argumente der Textur
    args = gvsParser->parse(args);

    // Durchsuche Parameter nach dem Typ der Texture
    std::string textureType;
    gvsParser->getParameter("type",textureType);

    if      ( textureType == "UniTex"       )      gvsP_init_unitex(gvsParser);
    else if ( textureType == "CheckerT2D"   )      gvsP_init_checkerT2D(gvsParser);
    //else if ( textureType == "ChequeredT2D" )      gvsP_init_chequeredT2D(gvsParser);
    //else if ( textureType == "ChequeredTimeT2D" )  gvsP_init_chequeredTimeT2D(gvsParser);
    else if ( textureType == "Image"        )      gvsP_init_image2dsampler(gvsParser);

    delete gvsParser;

    pointer R = ((sc->vptr->mk_symbol)(sc, "gtTexture"));
    return R;
}

/**
 * @brief gvsP_init_unitex
 * @param gP
 */
void gvsP_init_unitex (GvsParseScheme* gP) {

    GvsUniTex* unitex;
    double     col[3] = { 0.0, 0.0, 0.0 };

    if (gP->getParameter("colorgray",col)) {
        unitex = new GvsUniTex(col[0]);
    }

    if (gP->getParameter("color",col)) {
        unitex = new GvsUniTex(GvsColor(col[0],col[1],col[2]));
    }
    gpTexture.push_back(unitex);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-texture: ID already assigned!");
    }

    GvsTypeID tid = {gtTexture,gpTexture.size()-1,gpTexture[gpTexture.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
}

/**
 * @brief gvsP_init_checkerT2D
 * @param gP
 */
void gvsP_init_checkerT2D (GvsParseScheme* gP) {
    GvsCheckerT2D* texT2D;
    GvsTexture* currTex[2];

    gP->setAllowedName("texture",gp_string_string,0);

    int num = 0;
    int count = 0;
    bool texFound = false;
    std::string textureID;
    std::string msg;

    do {
        texFound = gP->getParameter("texture",textureID,num);
        if (texFound) {
            if (textureID=="gtTexture") {
                currTex[count] = gpTexture[gpTexture.size()-1];
            }
            else {
                if (gpTexture.size()>=1) {
                    getIDptr(gP,"CheckerT2D","Texture","texture",gtTexture,num);
                    currTex[count] = gpTexture[(gpTypeIDptr->second).vectorID];
                }
                else
                    scheme_error("init-shader: no texture available!");
            }
        }
        else {
            if (gpTexture.empty()) scheme_error("init-shader: no texture available!");
            else if (gpTexture.size()>1)
            {
                scheme_error("init-shader: texture-ID is missing!\n");
            }
            else
                currTex[count] = gpTexture[0];
#ifdef GVS_VERBOSE
            printf("init-shader: verwende eingegebene Textur\n");
#endif
        }
        //    currTex[count]->print(cerr);
        count++;
        num++;
    }
    while ( (count<2) );

    texT2D = new GvsCheckerT2D(currTex[0],currTex[1]);

    // Transformation
    m4d::Matrix<double,2,3> texMat2D; // = ScaleMat2D(2.0,2.0);
    if (!gP->getParameter("transform",&texMat2D)) texMat2D.setIdent();
    texT2D->setTransformation(texMat2D);

    gpTexture.push_back(texT2D);
     //texT2D->Print();

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-texture: ID already assigned!");
    }

    GvsTypeID tid = {gtTexture,gpTexture.size()-1,gpTexture[gpTexture.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));

}

//----------------------------------------------------------------------------
//         gvsP_init_chequeredT2D
//----------------------------------------------------------------------------
/*
void gvsP_init_chequeredT2D (GvsParseScheme* gP)
{
  GvsChequeredT2D* texT2D;
  GvsTexture* currTex[2];
  gP->setAllowedName("texture",gp_string_string,0);
  gP->setAllowedName("width",gp_string_double,1);

  int num = 0;
  int count = 0;
  bool texFound = false;
  string textureID;
  string msg;

  do
  {
    texFound = gP->getParameter("texture",textureID,num);

    if (texFound)
    {
      if (textureID=="gtTexture")
      {
        currTex[count] = gpTexture[gpTexture.size()-1];
      }
      else
      {
        if (gpTexture.size()>=1)
        {
          getIDptr(gP,"CheckerT2D","Texture","texture",gtTexture,num);
          currTex[count] = gpTexture[(gpTypeIDptr->second).vectorID];
        }
        else
          scheme_error("init-shader: keine Textur vorhanden!");
      }
    }
    else if (gpTexture.size()>1)
    {
      getIDptr(gP,"ChequeredT2D","Texture","texture",gtTexture,num);
      currTex[count] = gpTexture[(gpTypeIDptr->second).vectorID];
    }
    else
    {
      if (gpTexture.empty()) scheme_error("init-shader: keine Textur vorhanden!");
#ifdef GVS_VERBOSE
      printf("init-shader: verwende eingegebene Textur\n");
#endif
      currTex[count] = gpTexture[0];
    }
    //    currTex[count]->print(cerr);
    count++;
    num++;
  }
  while ( (count<2) );


  double width;
  if (!gP->getParameter("width",width)) width = 0.25;

  texT2D = new GvsChequeredT2D(currTex[0],currTex[1],width);

  // Transformation
  m4d::Matrix<double,2,3> texMat2D; // = ScaleMat2D(2.0,2.0);
  if (!gP->getParameter("transform",&texMat2D)) texMat2D.setIdent();
  texT2D->setTransformation(texMat2D);

  gpTexture.push_back(texT2D);

  string idname = "unknown";
  if (!gP->getParameter("id",idname))
  {
    appendNum(idname,gpTypeID.size());
  }
  else if (gpTypeID.find(idname)!=gpTypeID.end())
  {
    scheme_error("init-texture: ID schon vergeben!");
  }
  gpTypeID.insert(pair<string,GvsTypeID>(idname,(GvsTypeID)
  {
    gtTexture,gpTexture.size()-1,gpTexture[gpTexture.size()-1]
  }
                                        ));
}


//----------------------------------------------------------------------------
//         gvsP_init_chequeredT2D
//----------------------------------------------------------------------------
void gvsP_init_chequeredTimeT2D (GvsParseScheme* gP)
{
  GvsChequeredTimeT2D* texTimeT2D;
  GvsTexture* currTex[2];

  gP->setAllowedName("texture",gp_string_string,0);
  gP->setAllowedName("projector",gp_string_string,0);
  gP->setAllowedName("deltatbehaviour", gp_string_string,0);
  gP->setAllowedName("colorbehaviour",  gp_string_string,0);
  gP->setAllowedName("timebehaviour",   gp_string_string,0);
  gP->setAllowedName("deltat",gp_string_double,1);
  gP->setAllowedName("stretch", gp_string_double,1);


  int num = 0;
  int count = 0;
  bool texFound = false;
  string textureID;
  string msg;

  do
  {
    texFound = gP->getParameter("texture",textureID,num);

    if (texFound)
    {
      if (textureID=="gtTexture")
      {
        currTex[count] = gpTexture[gpTexture.size()-1];
      }
      else
      {
        if (gpTexture.size()>=1)
        {
          getIDptr(gP,"CheckerT2D","Texture","texture",gtTexture,num);
          currTex[count] = gpTexture[(gpTypeIDptr->second).vectorID];
        }
        else
          scheme_error("init-shader: keine Textur vorhanden!");
      }
    }
    else if (gpTexture.size()>1)
    {
      getIDptr(gP,"ChequeredT2D","Texture","texture",gtTexture,num);
      currTex[count] = gpTexture[(gpTypeIDptr->second).vectorID];
    }
    else
    {
      if (gpTexture.empty()) scheme_error("init-shader: keine Textur vorhanden!");
#ifdef GVS_VERBOSE
      printf("init-shader: verwende eingegebene Textur\n");
#endif
      currTex[count] = gpTexture[0];
    }
    //    currTex[count]->print(cerr);
    count++;
    num++;
  }
  while ( (count<2) );


  double deltat;
  if (!gP->getParameter("deltat",deltat)) deltat = 0.5;

  GvsProjector *proj = NULL;
  gpTypeIDptr = getIDptr(gP,"Device","Projector","projector",gtProjector);
  proj = gpProjector[(gpTypeIDptr->second).vectorID];


  texTimeT2D = new GvsChequeredTimeT2D(currTex[0],currTex[1],deltat, proj);

  string behaviour;
  if (gP->getParameter("deltatbehaviour",behaviour))
  {
    if (behaviour == string("deltat"))
      texTimeT2D->setBehaviour(useDeltaT);
    else if (behaviour == string("deltatlin"))
      texTimeT2D->setBehaviour(useDeltaTLin);
    else if (behaviour == string("scale"))
      texTimeT2D->setBehaviour(useScale);
    else
      scheme_error("init-texture: deltatbehaviour-Wert unbekannt!");
  }
  if (gP->getParameter("colorbehaviour",behaviour))
  {
    if (behaviour == string("specrend"))
      texTimeT2D->setBehaviour(useSpecrend);
    else if (behaviour == string("textures"))
      texTimeT2D->setBehaviour(useTextures);
    else
      scheme_error("init-texture: colorbehaviour-Wert unbekannt!");
  }
  if (gP->getParameter("timebehaviour",behaviour))
  {
    if (behaviour == string("coordtime"))
      texTimeT2D->setBehaviour(coordTime);
    else if (behaviour == string("localobstime"))
      texTimeT2D->setBehaviour(localObsTime);
    else if (behaviour == string("localintersectime"))
      texTimeT2D->setBehaviour(localIntersecTime);
    else
      scheme_error("init-texture: timebehaviour-Wert unbekannt!");
  }

  double stretch;
  if ( gP->getParameter("stretch",stretch) )
    texTimeT2D->setStretch(stretch);



  // Transformation
  m4d::Matrix<double,2,3> texMat2D; // = ScaleMat2D(2.0,2.0);
  if (!gP->getParameter("transform",&texMat2D)) texMat2D.setIdent();
  texTimeT2D->setTransformation(texMat2D);

  gpTexture.push_back(texTimeT2D);

  string idname = "unknown";
  if (!gP->getParameter("id",idname))
  {
    appendNum(idname,gpTypeID.size());
  }
  else if (gpTypeID.find(idname)!=gpTypeID.end())
  {
    scheme_error("init-texture: ID schon vergeben!");
  }
  gpTypeID.insert(pair<string,GvsTypeID>(idname,(GvsTypeID)
  {
    gtTexture,gpTexture.size()-1,gpTexture[gpTexture.size()-1]
  }
                                        ));
}
*/


//----------------------------------------------------------------------------
//         gvsP_init_image2dsampler
//----------------------------------------------------------------------------
void gvsP_init_image2dsampler ( GvsParseScheme* gP )
{
    GvsPicIOEnvelope* imgFile = new GvsPicIOEnvelope;
    GvsChannelImg2D*  chanImg = new GvsChannelImg2D;
    assert (chanImg!=NULL);

    gP->setAllowedName("file",gp_string_string,0);

    bool fileRead = false;
    std::string filename,msg;
    if (gP->getParameter("file",filename))
    {
        char file[256];
        strcpy(file,filename.c_str());
        fileRead = imgFile->readChannelImg(*chanImg,file);

        if (!fileRead) {
            msg = "Could not read file ";
            msg.append(filename);
            scheme_error(msg);
        }
    }

    GvsImg2DSampler* imgSampler = new GvsImg2DSampler(chanImg);

    m4d::Matrix<double,2,3> texMat2D;
    if (!gP->getParameter("transform",&texMat2D)) texMat2D.setIdent();
    imgSampler->setTransformation(texMat2D);


    gpTexture.push_back(imgSampler);

    std::string idname = "unknown";
    if (!gP->getParameter("id",idname)) {
        appendNum(idname,gpTypeID.size());
    }
    else if (gpTypeID.find(idname)!=gpTypeID.end()) {
        scheme_error("init-texture: ID already assigned!");
    }

    GvsTypeID tid = {gtTexture,gpTexture.size()-1,gpTexture[gpTexture.size()-1]};
    gpTypeID.insert(std::pair<std::string,GvsTypeID>(idname,tid));
    delete imgFile;
}
