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

#include "Shader/Surface/GvsSurfaceShader.h"
#include "Dev/GvsDevice.h"
#include "Img/GvsColor.h"
#include "Light/GvsLightSrc.h"
#include "Ray/GvsRay.h"
#include "Ray/GvsRayVisual.h"
#include "Texture/GvsTexture.h"


GvsSurfaceShader::GvsSurfaceShader() :
    intersection(NULL),
    matTexColor(NULL),
    matAmbient(1.0),
    matDiffuse(0.0) {
}

void GvsSurfaceShader::setPrimitiveColor( GvsTexture* tex ) {
    matTexColor = tex;
}

void GvsSurfaceShader::setKambient( double ka ) {
    matAmbient = ka;
}

void GvsSurfaceShader::setKdiffuse( double kl ) {
    matDiffuse = kl;
}


GvsColor GvsSurfaceShader::ambientColor ( GvsDevice* device, GvsRayVisual & ) {
    return device->lightSrcMgr->ambientLight();
}


GvsColor GvsSurfaceShader :: getIncidentLight( GvsDevice* device, GvsRayVisual& ray ) {
    if ( (ray.surfIntersec().getLocalDirection() | ray.surfIntersec().normNormal()) > 0.0 ) {
        ray.surfIntersec().invertNormal();
    }

    //GvsColor ambColor  = ambientColor( device, ray );
    GvsColor primColor = primitiveColor( ray.surfIntersec() );

    return primColor*matAmbient  +
            primColor*matDiffuse*reflArtifColor( device, ray ) +
            primColor*reflLightColor( device, ray );
}


GvsColor  GvsSurfaceShader::reflArtifColor( GvsDevice* device, GvsRayVisual &ray ) {    
    // normal and direction are given in standard object system
    m4d::vec3 intersecNormal = ray.surfIntersec().normNormal();
    m4d::vec3 inDirection    = -ray.surfIntersec().getLocalDirection();

    double cosNL = intersecNormal | inDirection;
    if (cosNL <= 0.0)
        return RgbBlack;

    GvsColor ambLight = device->lightSrcMgr->ambientLight();
    GvsColor outLight = cosNL * ambLight;
    return outLight;
}


GvsColor GvsSurfaceShader::reflLightColor( GvsDevice* device, GvsRayVisual &ray ) {
    GvsColor outLight = GvsColor(0.0);
    if (device->metric->getMetricName()!="Minkowski" ||
            !device->lightSrcMgr->withShadowRays()) {
        return outLight;
    }

    if (ray.surfIntersec().object()!=NULL) {
        if (ray.surfIntersec().object()->getObjType()!=inCoords ||
                ray.surfIntersec().object()->getMotion()!=NULL) {
            return outLight;
        }
    }

    // only valid for Minkowski spacetime in Cartesian coordinates !!!
    m4d::vec4 isecPoint = ray.surfIntersec().point();


    GvsLightSrc* light = NULL;
    bool lightFound = false;

    for(int i=0; i<(device->lightSrcMgr->length()); i++) {
        light = device->lightSrcMgr->getLightSrc(i,lightFound);

        if (lightFound) {
            m4d::vec4 lightPos = light->getPosition();
            m4d::vec3 lightDir = (lightPos.getAsV3D() - isecPoint.getAsV3D()).getNormalized();
            m4d::vec4 rayDir = m4d::vec4(-1.0,lightDir.x(0),lightDir.x(1),lightDir.x(2));

            // calculate a shadow ray
            GvsRayVisual* eyeRay = new GvsRayVisual(device->projector->getRayGen());
            eyeRay->recalc( isecPoint, rayDir );

            if (device->sceneGraph != NULL) {
                // find closest intersection of the shadow ray with a scene object
                if (!device->sceneGraph->testIntersection(*eyeRay)) {
                    outLight += light->color();
                }
            }

            delete eyeRay;
        }
    }
    return outLight;
}


GvsColor GvsSurfaceShader::primitiveColor( GvsSurfIntersec& intersec ) const {
    if (matTexColor != NULL) {
        return matTexColor->sampleColor( intersec );
    }
    return RgbBlack;
}


void GvsSurfaceShader::Print ( FILE* fptr ) {
    fprintf(fptr,"SurfaceShader {\n");
    fprintf(fptr,"\tobjcolor: "); matTexColor->Print(fptr);
    fprintf(fptr,"\tambient: %f\n",matAmbient);
    fprintf(fptr,"\tdiffuse: %f\n",matDiffuse);
    fprintf(fptr,"}\n");
}
