
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
//READING DONE
#include "GvsGramSchmidt.h"

#include <gsl/gsl_linalg.h>     //LINK:GSL
#include <gsl/gsl_matrix.h>		//LINK:GSL


GvsGramSchmidt :: GvsGramSchmidt( m4d::Metric* metric ) {
    mMetric = metric;
    for (int i=0; i<4; i++) {
        e[i] = m4d::vec4();
    }
    orthonormal = false;
}

GvsGramSchmidt :: GvsGramSchmidt ( m4d::Metric* metric, const m4d::vec4 pos) {
    mMetric = metric;
    setPosition(pos);
    for (int i=0; i<4; i++) {
        e[i] = m4d::vec4();
    }
    orthonormal = false;
}

GvsGramSchmidt :: GvsGramSchmidt ( m4d::Metric* metric, const m4d::vec4 pos,
                                   const m4d::vec4 e0, const m4d::vec4 u1, const m4d::vec4 u2, const m4d::vec4 u3 ) {
    mMetric = metric;
    setPosition(pos);
    setVectors (e0,u1,u2,u3);
    orthonormal = false;
}

GvsGramSchmidt :: ~GvsGramSchmidt() {
    mMetric = NULL;
}


void GvsGramSchmidt :: setMetric ( m4d::Metric* metric) {
    mMetric = metric;
}


void GvsGramSchmidt :: setPosition ( const m4d::vec4 pos ) {
    for (int i = 0; i < 5; i++)     {
        mPos[i] = pos.x(i);
    }
    assert( mMetric != NULL);
    mMetric->calculateMetric(mPos);
}


void GvsGramSchmidt :: setVectors  ( const m4d::vec4 e0, const m4d::vec4 u1, const m4d::vec4 u2, const m4d::vec4 u3 ) {
    e[0] = e0;
    normalizeVector(&e[0]);
    e[1] = u1;
    e[2] = u2;
    e[3] = u3;

    for (int i=0; i<4; i++) {
        if (e[i].isZero()) {
            std::cerr << "GramSchmidt:  Vektor " << i << " ist ungueltig!\n";
            exit(1);
        }
    }

    // test if orthonormal
    double sc;
    orthonormal = false;

    for (int i=0; i<4; i++) {			//PROBLEMATIC
        for (int j=i+1; j<4; j++) {
            sc = calcScalarProd ( &e[i], &e[j] );
            if (fabs(sc)<GVS_EPS)
                orthonormal |= true;
        }
    }
}


void GvsGramSchmidt :: getVectors  ( m4d::vec4 &e0, m4d::vec4 &e1, m4d::vec4 &e2, m4d::vec4 &e3 ) {
    e0 = e[0];
    e1 = e[1];
    e2 = e[2];
    e3 = e[3];
}

m4d::vec4 GvsGramSchmidt :: getVector ( int i ) const {
    return e[i];
}


bool GvsGramSchmidt :: isOrthonormal ( ) const {
    return orthonormal;
}


bool GvsGramSchmidt :: isRightHanded ( ) const {
    gsl_matrix* m = gsl_matrix_alloc(4,4);

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            gsl_matrix_set(m,i,j,e[i][j]);

    int signum;
    gsl_permutation *p = gsl_permutation_alloc(4);

    gsl_linalg_LU_decomp ( m, p, &signum );

    double det = gsl_linalg_LU_det(m,signum);

    gsl_permutation_free(p);
    gsl_matrix_free(m);

    return ((det>0) ? true : false);
}


void GvsGramSchmidt :: normalizeVector ( m4d::vec4* vec ) {
    double norm = fabs(calcScalarProd ( vec, vec ));
    //cerr << "GramSchmidt: Normalize Vektor:" << mMetric->getMetricName() << "-" << mMetric->getMetricCoeff(1,1) << "-" << mMetric->getMetricCoeff(2,2) << "-" << mMetric->getMetricCoeff(3,3) <<"\n";
    assert ( norm > 0.0 );

    norm = sqrt(norm);
    double val;
    for ( int i=0; i<4; i++ )     {
        val = (vec->x(i))/norm;
        if (fabs(val)<1.0e-15) val = 0.0;
        vec->setX(i,val);
    }
}


double GvsGramSchmidt :: calcScalarProd ( const m4d::vec4* v1, const m4d::vec4* v2 ) {
    double prod = 0.0;
    for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
            prod += mMetric->getMetricCoeff(i,j) * v1->x(i) * v2->x(j);
        }
    }
    return prod;
}


void GvsGramSchmidt :: calculateTetrad() {
    mMetric->calculateMetric(mPos);

    m4d::vec4 u1 = e[1];
    m4d::vec4 u2 = e[2];
    m4d::vec4 u3 = e[3];

    // projiziere auf e0

    e[1] = u1 + calcScalarProd(&e[0],&u1)*e[0];
    assert (!(e[1].isZero()));
    normalizeVector(&e[1]);

    e[2] = u2 + calcScalarProd(&e[0],&u2)*e[0] - calcScalarProd(&e[1],&u2)*e[1];
    assert (!(e[2].isZero()));
    normalizeVector(&e[2]);

    e[3] = u3 + calcScalarProd(&e[0],&u3)*e[0] - calcScalarProd(&e[1],&u3)*e[1] - calcScalarProd(&e[2],&u3)*e[2];
    assert (!(e[3].isZero()));
    normalizeVector(&e[3]);

    orthonormal = true;
}


void GvsGramSchmidt :: calculateTetrad ( m4d::vec4 &e0, m4d::vec4 &e1, m4d::vec4 &e2, m4d::vec4 &e3 ) {
    setVectors(e0,e1,e2,e3);
    calculateTetrad();
    getVectors(e0,e1,e2,e3);
}


void GvsGramSchmidt :: Print( FILE* fptr ) {
    fprintf(fptr,"GvsGramSchmidt {\n");
    fprintf(fptr,"\tpos:  %6.3f %6.3f %6.3f %6.3f\n",mPos[0],mPos[1],mPos[2],mPos[3]);
    fprintf(fptr,"\te0:   %6.3f %6.3f %6.3f %6.3f\n",e[0][0],e[0][1],e[0][2],e[0][3]);
    fprintf(fptr,"}\n");
}


void GvsGramSchmidt :: printS () const {
    printf("GvsGramSchmidt {\n");
    printf("\tPosition : ( %9.6f %9.6f %9.6f %9.6f %i )\n",mPos[0],mPos[1],mPos[2],mPos[3],0);
    printf("\te0 : ( %9.6f %9.6f %9.6f %9.6f )\n", e[0].x(0),e[0].x(1),e[0].x(2),e[0].x(3));
    printf("\te1 : ( %9.6f %9.6f %9.6f %9.6f )\n", e[1].x(0),e[1].x(1),e[1].x(2),e[1].x(3));
    printf("\te2 : ( %9.6f %9.6f %9.6f %9.6f )\n", e[2].x(0),e[2].x(1),e[2].x(2),e[2].x(3));
    printf("\te3 : ( %9.6f %9.6f %9.6f %9.6f )\n", e[3].x(0),e[3].x(1),e[3].x(2),e[3].x(3));
    printf("}\n");
    printf("rightRanded: %i\n",(int)isRightHanded());
}
