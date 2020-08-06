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

#include "Obj/GvsSceneObj.h"

GvsSceneObj::GvsSceneObj() {
    mObjType = local;
    mMetric  = NULL;
    mChart   = 0;
    stMotion  = NULL;
    haveMotion = false;
}

GvsSceneObj::GvsSceneObj(GvsObjType oTyp) {
    mObjType = oTyp;
    mMetric  = NULL;
    mChart   = 0;

    stMotion  = NULL;
    haveMotion = false;
}

GvsSceneObj::~GvsSceneObj()
{
    mMetric = NULL;
    stMotion = NULL;
}

void GvsSceneObj :: setObjType(GvsObjType oTyp) {
    mObjType = oTyp;
}

GvsObjType GvsSceneObj :: getObjType() const {
    return mObjType;
}

void GvsSceneObj::setChart(int chart ) {
    mChart = chart;
}

int GvsSceneObj::getChart() const {
    return mChart;
}


void GvsSceneObj::setMetric ( m4d::Metric* spacetime ) {
    mMetric = spacetime;
}

m4d::Metric* GvsSceneObj::getMetric () const {
    return mMetric;
}


void GvsSceneObj :: setMotion ( GvsStMotion *motion ) {
    if (stMotion!=NULL) delete stMotion;
    stMotion = motion;
    haveMotion = true;
}

GvsStMotion* GvsSceneObj :: getMotion ( void ) const {
    return stMotion;
}

/*
//----------------------------------------------------------------------------
//         stBoundBox
//----------------------------------------------------------------------------
void
GvsSceneObj :: calcSTBoundBox ( GvsBoundBox box )
{
  m4d::vec3 low = box.lowBounds();
  m4d::vec3 upp = box.uppBounds();

  // Bestimme die gr��ten Koordinatenwerte der normalen BoundingBox im lokalen System
  double m1 = GVS_MAX(low.x(0),upp.x(0));
  double m2 = GVS_MAX(low.x(1),upp.x(1));
  double m3 = GVS_MAX(low.x(2),upp.x(2));

  // Transformiere den Ortsvektor ins Koordinatensystem
  m4d::vec4 maxBound = m4d::vec4(1.0,m1,m2,m3);
  maxBound.print();
  m4d::vec4 maxBoundCoords = localTetrad->transToCoords(maxBound);

  // Bestimme daraus den maximalen raumartigen und zeitartigen Abstand der BoundingBox
  GvsMetric* metric = localTetrad->getMetric();
  double spaceDist,timeDist;
  metric->calcDist(localTetrad->getPosition(),maxBoundCoords,spaceDist,timeDist);

  cout << "Distances: " << spaceDist << " " << timeDist << std::endl;

  // Setze hier spaceDist und timeDist als stBoundingBox
  //  stBoundBox = GvsBoundBox4D(minBoundCoords,maxBoundCoords);
}
*/


bool GvsSceneObj :: testIntersection(GvsRay&) {
    std::cerr << "Error in GvsSceneObj::testIntersection(GvsRay&): not implemented." << std::endl;
    return false;
}

/*
bool
GvsSceneObj::testLocalIntersection( GvsRay &ray, const int rayPartIndex, const int seg,
                                    GvsLocalTetrad* lt0, GvsLocalTetrad* lt1,
                                    const m4d::vec4 p0, const m4d::vec4 p1,
                                    double& tau0, double& tau1 )
{
  std::cerr << "Error in GvsSceneObj::testLocalIntersection(GvsRay&): not implemented." << std::endl;
  return false;
}
*/

bool
GvsSceneObj::testLocalIntersection(GvsRay &, const int ,
                                    GvsLocalTetrad* , GvsLocalTetrad* ,
                                    const m4d::vec4 , const m4d::vec4 ,
                                    double& , double&  )
{
    std::cerr << "Error in GvsSceneObj::testLocalIntersection(GvsRay&): not implemented." << std::endl;
    return false;
}
