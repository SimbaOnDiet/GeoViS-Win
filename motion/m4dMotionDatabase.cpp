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
// -------------------------------------------------------------------------------
/*
   m4dMetricDatabase.cpp

  Copyright (c) 2009-2014-2011  Thomas Mueller, Frank Grave


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/
// -------------------------------------------------------------------------------

#include "m4dMotionDatabase.h"

namespace m4d
{

/*!
 */
IntegratorDatabase :: IntegratorDatabase( bool printDatabase )
{
    init();
    if (printDatabase)
        printIntegratorList();
}

IntegratorDatabase :: ~IntegratorDatabase()
{
    if (!mIntegratorMap.empty())
        mIntegratorMap.clear();
    if (!mIntegratorNickMap.empty())
        mIntegratorNickMap.clear();
}


// *********************************** public methods ******************************
/*! Get the number of implemented metrics.
 *
 *  \return int: number of implemented metrics.
 */
int 
IntegratorDatabase :: getNumIntegrators ( )
{
    return NUM_GEOD_SOLVERS;
}

/*!  Initialize metric 'num' and return the pointer to it.
 *
 *  \param cMetric : pointer to metric.
 *  \param num     : metric number.
 */
Geodesic*
IntegratorDatabase :: getIntegrator (Metric* cMetric, enum_integrator num )
{
    return initializeIntegrator(cMetric, num);
}

/*! Initialize metric 'mName' and return the pointer to it.
 *
 *  \param cMetric : pointer to metric.
 *  \param mName   : name of metric.
 */
Geodesic*
IntegratorDatabase :: getIntegrator (Metric* cMetric, std::string mName )
{
    mIntegratorMapItr = mIntegratorNickMap.find( mName );
    if (mIntegratorMapItr == mIntegratorMap.end())
    {
        mIntegratorMapItr = mIntegratorMap.find ( mName );
        if (mIntegratorMapItr == mIntegratorMap.end()) {
            fprintf(stderr,"Integrator '%s' is not implemented!\n", mName.c_str());
            return NULL;
        }
    }

    return initializeIntegrator(cMetric, mIntegratorMapItr->second );
}

/*! Get the name of metric 'num'.
 *
 *  \param num : number of metric.
 *  \return string : name of metric.
 */
std::string 
IntegratorDatabase :: getIntegratorName ( enum_integrator num )
{
    if ( int(num) >= 0 && int(num) < NUM_GEOD_SOLVERS )
        return stl_solver_names[num];

    return std::string();
}

/*! Get the number of the 'mName' metric.
 *
 *  \param mName : name of metric.
 *  \return enum_metric : number of metric.
 */
enum_integrator
IntegratorDatabase :: getIntegratorNr( std::string mName )
{
    mIntegratorMapItr = mIntegratorNickMap.find( mName );
    if (mIntegratorMapItr == mIntegratorNickMap.end())
    {
        mIntegratorMapItr = mIntegratorMap.find ( mName );
        if (mIntegratorMapItr == mIntegratorMap.end()) {
            fprintf(stderr,"Integrator '%s' is not implemented!\n", mName.c_str());
            return gsUnknown;
        }
    }

    return mIntegratorMapItr->second;
}

/*!  Print list of all metrics.
 *
 *  \param fptr : pointer to file.
 */
void 
IntegratorDatabase :: printIntegratorList ( FILE* fptr )
{
    fprintf(fptr,"  Solver nicknanme            Solver full name\n");
    fprintf(fptr,"-----------------------------------------------------------\n");
    for(int i=0; i<NUM_GEOD_SOLVERS; i++) {
        fprintf(fptr,"%-23s : %s\n",stl_solver_nicknames[i],stl_solver_names[i]);
    }
}


// ******************************** protected methods ******************************

/*!  Initialize database: store the metrics in a map.
 *
 */
void 
IntegratorDatabase :: init ( )
{
    for(int i=0; i<NUM_GEOD_SOLVERS; i++) {
        mIntegratorMap.insert( std::pair<std::string,enum_integrator>(stl_solver_names[i], enum_integrator(i)) );
        mIntegratorNickMap.insert( std::pair<std::string,enum_integrator>(stl_solver_nicknames[i], enum_integrator(i)) );
    }
}


// ---------------------------------------------------
//           Please sort by enum name !!
// ---------------------------------------------------
/*! Initialize metric: generate a new instance.
 *
 *  \param cMetric : pointer to metric.
 *  \param num     : enumeration of metric.
 */
Geodesic*
IntegratorDatabase :: initializeIntegrator (Metric* cMetric,  enum_integrator  num )
{

    Geodesic*  currGeo = NULL;

    const gsl_odeiv_step_type*  step_type = gsl_odeiv_step_rk4;
    switch (num)
    {
        case m4d::gsIrk4:
#ifdef USE_DP_INT
        case m4d::gsIdp54:
        case m4d::gsIdp65:
#endif
        case gsInrbs:
            break;
        case m4d::gsIgslrk2:
            step_type = gsl_odeiv_step_rk2;
            break;
        case m4d::gsIgslrk4:
            step_type = gsl_odeiv_step_rk4;
            break;
        case m4d::gsIgslfehlberg:
            step_type = gsl_odeiv_step_rkf45;
            break;
        case m4d::gsIgslcash:
            step_type = gsl_odeiv_step_rkck;
            break;
        case m4d::gsIgslprinc:
            step_type = gsl_odeiv_step_rk8pd;
            break;
        case m4d::gsIi2:
            step_type = gsl_odeiv_step_rk2imp;
            break;
            /*
    case m4d::gsIbs:
      step_type = gsl_odeiv_step_bsimp;
      break;
      */
        case m4d::gsIm1:
            step_type = gsl_odeiv_step_gear1;
            break;
        case m4d::gsIm2:
            step_type = gsl_odeiv_step_gear2;
            break;
        default: break;
    }

	//bugfix_step_type = step_type;

    if (num==gsIrk4)
        currGeo = new GeodesicRK4( cMetric);
    else if (num==gsInrbs) {
        currGeo = new GeodesicBS( cMetric );
    }
#ifdef USE_DP_INT
    else if (num==gsIdp54) {
        currGeo = new GeodesicDP54( cMetric );
    }
    else if (num==gsIdp65) {
        currGeo = new GeodesicDP65( cMetric );
    }
#endif
    else {
        currGeo = new GeodesicGSL( cMetric, step_type, num);
    }

    // mData->geodSolver->setBoundingBox( m4d::vec4(-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX), m4d::vec4(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX) );
    // mData->geodSolver->setGeodesicType( mData->type );

    return currGeo;
}

} // end namespace m4d
