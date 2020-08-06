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
#ifndef GVS_PARSER_H
#define GVS_PARSER_H

extern "C" {
#include "scheme-private.h"
#include "scheme.h"
}

#include <Parser/parse_camera.h>
#include <Parser/parse_compObj.h>
#include <Parser/parse_device.h>
#include <Parser/parse_helper.h>
#include <Parser/parse_locTed.h>
#include <Parser/parse_light.h>
#include <Parser/parse_metric.h>
#include <Parser/parse_motion.h>
#include <Parser/parse_planeObjects.h>
#include <Parser/parse_projector.h>
#include <Parser/parse_raygen.h>
#include <Parser/parse_shader.h>
#include <Parser/parse_solidObjects.h>
#include <Parser/parse_csgObj.h>

#include <Parser/parse_solver.h>
#include <Parser/parse_texture.h>
#include <Parser/parse_transformation.h>
#include <Parser/parse_transformObject.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>

#include "GvsGlobalDefs.h"
#include "Cam/GvsCamera.h"
#include "Dev/GvsDevice.h"
#include "Dev/GvsProjector.h"
#include "Img/GvsChannelImg2D.h"
#include "Light/GvsLightSrc.h"
#include "Light/GvsLightSrcMgr.h"
#include "Obj/STMotion/GvsLocalTetrad.h"
#include "Obj/GvsSceneObj.h"
#include "Obj/Comp/GvsLocalCompObj.h"
#include "Obj/STMotion/GvsStMotion.h"
#include "Obj/SolidObj/GvsSolidCSGObj.h"
#include "Obj/SolidObj/GvsSolidDifferObj.h"
#include "Obj/SolidObj/GvsSolidIntersecObj.h"
#include "Obj/SolidObj/GvsSolidUnifiedObj.h"
#include "Ray/GvsRayGen.h"
#include "Shader/GvsShader.h"
#include "Utils/GvsGeodSolver.h"

//#include <Texture/GvsTexture.h>

#include <metric/m4dMetricDatabase.h>

class GvsDevice;

/**
 * @brief The GvsParser class
 */
class GvsParser
{
public:
    GvsParser();
    virtual ~GvsParser();

    void  deleteAll();
    void  read_scene     ( const char* name );
    int   getNumDevices  ( void ) const;

    void  initStandard   ( GvsDevice *device );
    void  getDevice      ( GvsDevice *device, unsigned int k = 0 );

    void  printAll       ( FILE *fptr = stderr ) const;


    m4d::Metric*      getMetric    ( unsigned int k = 0 );
    GvsGeodSolver*    getSolver    ( unsigned int k = 0 );
    GvsCamera*        getCamera    ( unsigned int k = 0 );
    GvsProjector*     getProjector ( unsigned int k = 0 );
    GvsRayGen*        getRayGen    ( unsigned int k = 0 );
    GvsLightSrcMgr*   getLightMgr  ( unsigned int k = 0 );
    GvsLightSrc*      getLight     ( unsigned int k = 0 );
    GvsSceneObj*      getSceneObj  ( unsigned int k = 0 );
    GvsStMotion*      getMotion    ( unsigned int k = 0 );

protected:
    std::string  getFullPathname();
};


#endif
