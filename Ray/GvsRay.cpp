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

#include "Ray/GvsRay.h"
#include "Ray/GvsRayGen.h"


GvsRay :: GvsRay() {
    rayGen   = NULL;
    rayPoints= NULL;
    rayDirs  = NULL;
    rayLambda = NULL;
    raySachs1 = NULL;
    raySachs2 = NULL;
    rayJacobi = NULL;

    rayTetrad = NULL;
    rayHasTetrad = false;

    rayID = getNextRayID();

    rayNumPoints = 0;
}

GvsRay :: GvsRay ( GvsRayGen* gen ) {
    assert (gen!=NULL);

    rayGen = gen;
    rayPoints = NULL;
    rayDirs   = NULL;
    rayLambda = NULL;
    raySachs1 = NULL;
    raySachs2 = NULL;
    rayJacobi = NULL;

    rayTetrad = NULL;
    rayHasTetrad = false;

    rayID = getNextRayID();

    rayNumPoints = 0;
}

GvsRay :: GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, GvsRayGen* gen ) {
    assert (gen!=NULL);
    rayGen = gen;
    rayLambda = NULL;
    raySachs1 = NULL;
    raySachs2 = NULL;
    rayJacobi = NULL;

    rayID = getNextRayID();

    rayDirs = new m4d::vec4[rayGen->getMaxNumPoints()];

    rayNumPoints = 0;
    rayPoints= rayGen->calcPolyline(orig,dir,rayDirs,rayNumPoints);   // construct ray
    assert(rayNumPoints >= 2);

    rayMinSearchDist = GVS_EPS;
    rayMaxSearchDist = double(rayNumPoints-1); // ???


    rayTetrad = NULL;
    rayHasTetrad = false;
}


GvsRay :: GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, GvsRayGen* gen,
                   double minSearchDist, double maxSearchDist )
{
    // std::cerr << "GvsRay :: GvsRay()...min.max.\n";
    assert (gen!=NULL);
    rayGen = gen;
    rayLambda = NULL;
    raySachs1 = NULL;
    raySachs2 = NULL;
    rayJacobi = NULL;

    rayID = getNextRayID();

    rayDirs = new m4d::vec4[rayGen->getMaxNumPoints()];

    rayNumPoints = 0;
    rayPoints= rayGen->calcPolyline(orig,dir,rayDirs,rayNumPoints);

    setSearchInterval(minSearchDist,maxSearchDist);

    rayTetrad = NULL;
    rayHasTetrad = false;
}



GvsRay :: GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad, GvsRayGen* gen )
{
    // std::cerr << "GvsRay :: GvsRay()...\n";
    assert (gen!=NULL);
    rayGen = gen;
    rayLambda = NULL;
    raySachs1 = NULL;
    raySachs2 = NULL;
    rayJacobi = NULL;

    rayID = getNextRayID();

    rayPoints= NULL;
    rayDirs  = NULL;

    rayNumPoints = 0;
    rayTetrad = rayGen->calcParTransport(orig,dir,tetrad,rayNumPoints);
    assert(rayNumPoints >= 2);
    rayHasTetrad = true;

    rayMinSearchDist = GVS_EPS;
    rayMaxSearchDist = double(rayNumPoints-1); // ???
}


GvsRay :: GvsRay ( const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad, GvsRayGen* gen,
                   double minSearchDist, double maxSearchDist )
{
    assert (gen!=NULL);
    rayGen = gen;

    rayID = getNextRayID();

    rayPoints = NULL;
    rayDirs   = NULL;
    rayLambda = NULL;
    raySachs1 = NULL;
    raySachs2 = NULL;
    rayJacobi = NULL;

    rayNumPoints = 0;
    rayTetrad = rayGen->calcParTransport(orig,dir,tetrad,rayNumPoints);
    assert(rayNumPoints >= 2);
    rayHasTetrad = true;

    setSearchInterval(minSearchDist,maxSearchDist);
}


GvsRay::~GvsRay()
{
    rayGen = NULL;
    deleteAll();
    rayHasTetrad = false;
}

void GvsRay::deleteAll() {

    if (rayPoints!=NULL) {
        delete [] rayPoints;
        rayPoints = NULL;
    }

    if (rayDirs!=NULL) {
        delete [] rayDirs;
        rayDirs = NULL;
    }

    if (rayTetrad!=NULL) {
        delete [] rayTetrad;
        rayTetrad = NULL;
    }

    if (rayLambda!=NULL) {
        delete [] rayLambda;
        rayLambda = NULL;
    }

    if (raySachs1!=NULL) {
        delete [] raySachs1;
        raySachs1 = NULL;
    }

    if (raySachs2!=NULL) {
        delete [] raySachs2;
        raySachs2 = NULL;
    }

    if (rayJacobi!=NULL) {
        delete [] rayJacobi;
        rayJacobi = NULL;
    }
}


void GvsRay :: setSearchInterval ( double minDist, double maxDist ) {
    if ( minDist < 0.0 ) {
        minDist = 0.0;
    }
    if ( maxDist < 0.0 ) {
        maxDist = 0.0;
    }

    if ( minDist < maxDist ) {
        rayMinSearchDist = minDist;
        rayMaxSearchDist = maxDist;
    }
    else {
        rayMinSearchDist = maxDist;
        rayMaxSearchDist = minDist;
    }
}

void GvsRay :: setMinSearchDist ( double minDist ) {
    assert ( minDist >= 0.0 );
    rayMinSearchDist = minDist;

    if (rayMinSearchDist > rayMaxSearchDist)
    {
        double t = rayMaxSearchDist;
        rayMaxSearchDist = rayMinSearchDist;
        rayMinSearchDist = t;
    }
}

void
GvsRay :: setMaxSearchDist ( double maxDist )
{
    assert ( maxDist >= 0.0 );
    rayMaxSearchDist = maxDist;

    if (rayMinSearchDist > rayMaxSearchDist)
    {
        double t = rayMaxSearchDist;
        rayMaxSearchDist = rayMinSearchDist;
        rayMinSearchDist = t;
    }
}



bool GvsRay :: recalc ( const m4d::vec4 &orig, const m4d::vec4 &dir ) {
    assert (rayGen != NULL);

    deleteAll();
    rayHasTetrad = false;

    rayID = getNextRayID();

    rayNumPoints = 0;
    rayPoints = rayGen->calcPolyline(orig,dir,rayNumPoints);

    if (rayPoints==NULL || rayNumPoints<2) {
        return false;
    }

    rayMinSearchDist = GVS_EPS;
    rayMaxSearchDist = double(rayNumPoints-1);
    return true;
}



bool GvsRay :: recalc (const m4d::vec4 &orig, const m4d::vec4 &dir, const GvsLocalTetrad *tetrad ) {
    // std::cerr << "GvsRay::recalc()...\n";
    //  tetrad.print(cerr);

    assert (rayGen != NULL);
    deleteAll();

    rayNumPoints = 0;
    //rayGen->Print();

    rayTetrad = rayGen->calcParTransport(orig,dir,tetrad,rayNumPoints);
    if (rayTetrad==NULL || rayNumPoints<2) {
        return false;
    }
    rayHasTetrad = true;

    rayMinSearchDist = GVS_EPS;
    rayMaxSearchDist = double(rayNumPoints-1);
    return true;
}


bool GvsRay::recalcJacobi ( const m4d::vec4 &orig, const m4d::vec4 &dir,
                            const m4d::vec3 &locRayDir, const GvsLocalTetrad* tetrad ) {
    assert (rayGen != NULL);
    //fprintf(stderr,"GvsRay::recalcJacobi() ... \n");
    deleteAll();

    rayNumPoints = 0;
    rayGen->calcSachsJacobi(orig,dir,locRayDir,tetrad,
                            rayPoints,rayDirs,rayLambda,raySachs1,raySachs2,rayJacobi,rayMaxJacobi,rayNumPoints);
    if (rayNumPoints<2) {
        return false;
    }

    rayMinSearchDist = GVS_EPS;
    rayMaxSearchDist = double(rayNumPoints-1);
    return true;
}



GvsRayGen*  GvsRay::getRayGen() const {
    return rayGen;
}


m4d::vec4* GvsRay::points() {
    assert(rayNumPoints>=2);

    if ((rayHasTetrad) && (rayPoints==NULL)) {
        rayPoints= new m4d::vec4[rayNumPoints];
        for (int i=0; i<rayNumPoints; i++)
            rayPoints[i] = rayTetrad[i].getPosition();
    }
    return rayPoints;
}


m4d::vec4* GvsRay::tangents() {
    return rayDirs;
}

int GvsRay :: getNumPoints() const {
    return rayNumPoints;
}

m4d::vec4 GvsRay :: getPoint ( int index ) const {
    assert ( (index >= 0) && (index < rayNumPoints) );
    if (rayHasTetrad && rayTetrad!=NULL) {
        return rayTetrad[index].getPosition();
    } else {
        if (rayPoints!=NULL) {
            return rayPoints[index];
        }
    }
    return m4d::vec4(0);
}

m4d::vec4 GvsRay :: getTangente ( int index ) const {
    assert ( (index >= 0) && (index < rayNumPoints) );
    if (rayHasTetrad && rayTetrad!=NULL) {
        return rayTetrad[index].getVelocity();
    } else {
        if (rayDirs!=NULL) {
            return rayDirs[index];
        }
    }
    return m4d::vec4(0);
}

m4d::vec5 GvsRay::getJacobi( int index ) const {
    assert ( (index >= 0) && (index < rayNumPoints) );
    if (rayJacobi!=NULL) {
        return rayJacobi[index];
    }
    return m4d::vec5(0);
}

GvsLocalTetrad GvsRay :: getTetrad ( int index ) const {
    assert ( (index >= 0) && (index < rayNumPoints) );
    if (rayHasTetrad && rayTetrad!=NULL) {
        return rayTetrad[index];
    }
    return GvsLocalTetrad();
}

void GvsRay :: setPoints( m4d::vec4* points ) {
    if ( rayPoints!= NULL ) {
        delete [] rayPoints;
    }
    rayPoints = points;
}

void GvsRay :: setDirs( m4d::vec4* dirs ) {
    if ( rayDirs != NULL ) {
        delete [] rayDirs;
    }
    rayDirs = dirs;
}

void GvsRay :: setNumPoints( int noPts ) {
    rayNumPoints = noPts;

    rayMinSearchDist = GVS_EPS;
    rayMaxSearchDist = double(rayNumPoints-1); // ???
}

ulong GvsRay :: getID () const {
    return rayID;
}


double GvsRay :: minSearchDist ( ) const {
    return rayMinSearchDist;
}

double GvsRay :: maxSearchDist ( ) const {
    return rayMaxSearchDist;
}

bool GvsRay :: isValidSurfIntersec ( double dist ) const {
    return (( dist > rayMinSearchDist ) && ( dist < rayMaxSearchDist ));
}

ulong GvsRay :: getNextRayID() {
    static ulong RayID = 1UL;

    if ( RayID == ULONG_MAX )
        return (RayID = 1UL);
    else
        return RayID++;
} 

void GvsRay :: timeShiftRay( double timeDelta ) {
    if (!rayHasTetrad) {
        assert(rayPoints!=NULL);
        for (int i=0; i<rayNumPoints; i++)
            rayPoints[i].setX( 0, rayPoints[i].x(0)+timeDelta);
    }
    else {
        assert(rayTetrad!=NULL);
        for (int i=0; i<rayNumPoints; i++)
            rayTetrad[i].setPositionX( 0, rayTetrad[i].getPosition().x(0)+timeDelta );
    }
}


void GvsRay :: Print( FILE* fptr ) {
    if (!rayHasTetrad)     {
        assert(rayPoints!=NULL);
        for (int i=0; i<rayNumPoints; i++) {
            rayPoints[i].printS(fptr);
        }
    }
    else {
        assert(rayTetrad!=NULL);

        m4d::vec4 pos;
        m4d::vec4 vel;
        m4d::vec4 base[4];

        for (int i=0; i<rayNumPoints; i++)
        {
            pos = rayTetrad[i].getPosition();
            vel = rayTetrad[i].getVelocity();
            for (int j=0; j<4; j++)
                base[j] = rayTetrad[i].getE(j);

            for (int k=0; k<4; k++) {
                fprintf(fptr,"%10.6f ",pos[k]);
            }
            fprintf(fptr," ");
            for (int k=0; k<4; k++) {
                fprintf(fptr,"%10.6f ",vel[k]);
            }
            fprintf(fptr," ");
            for (int k=0; k<4; k++) {
                for (int j=0; j<4; j++) {
                    fprintf(fptr,"%10.6f ",base[k][j]);
                }
            }
            fprintf(fptr,"\n");
        }
    }
}


double GvsRay::calcRayDist( const int seg, const double alpha ) {
    return double(seg)+alpha;
}

bool GvsRay::isIn( const int seg, const double alpha, const double endDist ) {
    return ( GvsRay::isValidAlpha(alpha) &&
             GvsRay::calcRayDist(seg,alpha) >= 0 &&
             GvsRay::calcRayDist(seg,alpha) <  endDist );
}

bool GvsRay::isValidAlpha( const double alpha ) {
    return (alpha >= 0.0) && (alpha < 1.0);
}

void GvsRay::splitRayDist( const double distance, int &seg, double &alpha ) {
    alpha = distance-floor(distance);
    seg = int(distance);
}
