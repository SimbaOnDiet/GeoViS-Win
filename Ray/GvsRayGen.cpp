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

#include "Ray/GvsRayGen.h"

GvsRayGen :: GvsRayGen() :
    actualSolver(NULL) {
    maxNumPoints = 2;
}

GvsRayGen :: GvsRayGen( GvsGeodSolver *solver ) :
    actualSolver(solver) {
    actualSolver->setGeodType(m4d::enum_geodesic_lightlike);
    actualSolver->setTimeDir(m4d::enum_time_backward);
    maxNumPoints = 2;
}

GvsRayGen :: GvsRayGen( GvsGeodSolver *solver,
                        m4d::enum_geodesic_type type, m4d::enum_time_direction dir) :
    actualSolver(solver)  {
    actualSolver->setGeodType(type);
    actualSolver->setTimeDir(dir);
    maxNumPoints = 2;
}

GvsRayGen :: ~GvsRayGen() {
    actualSolver = NULL;
}


void GvsRayGen :: setBoundBox ( const GvsBoundBox4D &box ) {
    boundBox = box;
    actualSolver->setBoundingBox(box.lowBounds(),box.uppBounds());
}

GvsBoundBox4D GvsRayGen :: getBoundBox ( ) const {
    return boundBox;
}

void GvsRayGen :: setMaxNumPoints(const int maxPoints) {
    maxNumPoints = maxPoints;
}

int GvsRayGen :: getMaxNumPoints() const {
    return maxNumPoints;
}


void GvsRayGen :: setActualSolver( GvsGeodSolver *solver ) {
    actualSolver = solver;
    actualSolver->setGeodType(m4d::enum_geodesic_lightlike);
    actualSolver->setTimeDir(m4d::enum_time_backward);
}


GvsGeodSolver* GvsRayGen :: getActualSolver ( ) const {
    return actualSolver;
}

/**
 * @param startOrig   Initial position of the ray in coordinates.
 * @param startDir    Initial direction of the ray in coordinates.
 * @param numPoints   Reference to number of points calculated.
 * @return pointer to ray points.
 */
m4d::vec4* GvsRayGen :: calcPolyline(const m4d::vec4 &startOrig, const m4d::vec4 &startDir, int &numPoints) {    
    assert(actualSolver!=NULL);

    m4d::Metric* metric = actualSolver->getMetric();
    if ( metric!=NULL) {
        if (metric->breakCondition(startOrig) ) {
            std::cerr << "Error in GvsRayGenSimple :: calcPolyline" << std::endl;
            std::cerr << "StartPos already satisfies breakCondition" << std::endl;
            return NULL;
        }
    }

    m4d::vec4* points = NULL;
    m4d::vec4* dirs = NULL;
    actualSolver->calculateGeodesic(startOrig,startDir,maxNumPoints,points,dirs,numPoints);
    delete [] dirs;
    return points;
}


m4d::vec4* GvsRayGen :: calcPolyline(const m4d::vec4 &startOrig, const m4d::vec4 &startDir,
                                     m4d::vec4 *&dirs, int &numPoints)
{
    // startOrig and startDir have to be given in coordinates !!    
    assert(actualSolver!=NULL);

    m4d::Metric* metric = actualSolver->getMetric();
    if ( metric!=NULL) {
        if (metric->breakCondition(startOrig) ) {
            std::cerr << "error in GvsRayGenSimple :: calcPolyline" << std::endl;
            std::cerr << "StartPos already satisfies breakCondition" << std::endl;
            return NULL;
        }
    }

    m4d::vec4* points = NULL;
    actualSolver->calculateGeodesic(startOrig,startDir,maxNumPoints,points,dirs,numPoints);
    return points;
}


//----------------------------------------------------------------------------
//         calcParTransport
//----------------------------------------------------------------------------
//  Neben der Geodaetenintegration wird noch zusaetzlich die lokale Tetrade
//  des Beobachters bis zum Emissionszeitpunkt parallel-transportiert.
//
GvsLocalTetrad*
GvsRayGen :: calcParTransport ( const m4d::vec4 &startOrig, const m4d::vec4 &startDir, const GvsLocalTetrad *lt,
                                int &numPoints )
{
    //  std::cerr << "GvsRayGenSimple :: calcParTransport () ...\n";
    assert(actualSolver!=NULL);

    m4d::Metric* metric = actualSolver->getMetric();
    if ( metric!=NULL) {
        if (metric->breakCondition(startOrig) ) {
            std::cerr << "error in GvsRayGenSimple :: calcPolyline" << std::endl;
            std::cerr << "StartPos already satisfies breakCondition" << std::endl;
            return NULL;
        }
    }

    m4d::vec4 base[4]; // = new m4d::vec4[4];
    for (int i=0; i<4; i++) {
        base[i] = lt->getE(i);
    }

    GvsLocalTetrad* localT = NULL;
    actualSolver->calcParTransport(startOrig,startDir,base,maxNumPoints,localT,numPoints);
    return localT;
}


bool  GvsRayGen::calcSachsJacobi ( const m4d::vec4 &startOrig, const m4d::vec4 &startDir, const m4d::vec3 &localDir,
                                   const GvsLocalTetrad* lt,
                                   m4d::vec4 *&points, m4d::vec4 *&dirs, double *&lambda,
                                   m4d::vec4 *&sachs1, m4d::vec4 *&sachs2,
                                   m4d::vec5 *&rayJacobi, m4d::vec5 &rayMaxJacobi,
                                   int &numPoints ) {
    assert(actualSolver!=NULL);

    m4d::Metric* metric = actualSolver->getMetric();
    if ( metric!=NULL) {
        if (metric->breakCondition(startOrig) ) {
            std::cerr << "Error in GvsRayGenSimple :: calcPolyline" << std::endl;
            std::cerr << "StartPos already satisfies breakCondition" << std::endl;
            return false;
        }
    }
    actualSolver->calcSachsJacobi(startOrig,startDir,localDir,maxNumPoints,
                                  lt,points,dirs,lambda,sachs1,sachs2,rayJacobi,rayMaxJacobi,numPoints);
    return true;
}



void GvsRayGen::Print( FILE* fptr )
{
    fprintf(fptr,"RayGen {\n");
    fprintf(fptr,"\tmaxNumPoints: %d",maxNumPoints);
    fprintf(fptr,"}\n");
}
