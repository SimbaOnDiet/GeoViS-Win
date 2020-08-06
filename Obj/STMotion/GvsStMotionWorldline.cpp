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
//  Copyright (c) 2013, Universitaet Stuttgart, VISUS, Thomas Mueller
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

#include "GvsStMotionWorldline.h"
#include <fstream>

using namespace std;

//----------------------------------------------------------------------------
//         constructor, destructor
//----------------------------------------------------------------------------
GvsStMotionWorldline::GvsStMotionWorldline()
    : GvsStMotion()
{
  mSpacetime = NULL;
  mSolver    = NULL;
  mType      = gvsMotionGeneral;

  lambdaStep = 0.01;
  maxNumPoints = 10000;
}

GvsStMotionWorldline::GvsStMotionWorldline(m4d::Metric* spacetime)
    : GvsStMotion()
{
//    cerr << "Initialize GvsStMotionWorldline..." << endl;
  mSpacetime = spacetime;
  mSolver    = NULL;
  mType      = gvsMotionGeneral;

  lambdaStep = 0.01;
  maxNumPoints = 10000;

  if (!localTetrad.empty()) deleteAllEntries();
}

GvsStMotionWorldline::GvsStMotionWorldline(m4d::Metric* spacetime,  GvsGeodSolverBase* solver)
    : GvsStMotion()
{
//    cerr << "Initialize GvsStMotionWorldline..." << endl;
  mSpacetime = spacetime;
  mSolver    = solver;
  mType      = gvsMotionGeneral;

  lambdaStep = 0.01;
  maxNumPoints = 10000;

  if (!localTetrad.empty()) deleteAllEntries();

  mSolver->setGeodCharacter(timelike,forward,gvsDestCoDir);
}

GvsStMotionWorldline::~GvsStMotionWorldline()
{
}

//----------------------------------------------------------------------------
//         setStartLocalTetrad
//----------------------------------------------------------------------------
void
GvsStMotionWorldline :: setStartLT ( GvsLocalTetrad *lt )
{
  localTetrad.push_back(lt);
  numPositions = localTetrad.size();
}

//----------------------------------------------------------------------------
//         setStartVal
//----------------------------------------------------------------------------
void
GvsStMotionWorldline :: setStartVal ( const P5D &pos, V4D &vel, bool inCoords, GvsLFType lfType )
{
  assert (numPositions==0);

  V4D coordVel = vel;

  // Ist die Geschwindigkeit nicht in Koordinaten gegeben, so muss diese noch umgewandelt werden.
  // Die Transformation ist abhaengig davon, bzgl welcher Tetrade die Geschw. gegeben ist.
  if (!inCoords)
  {
    // Geschwindigkeit ist bzgl lokaler Tetrade gegeben.

    double vq = vel.x(1)*vel.x(1) + vel.x(2)*vel.x(2) + vel.x(3)*vel.x(3);
    assert ( vq < 1.0 );// Pruefe ob v<c
    double gamma = 1.0/sqrt(1.0-vq);

    // lokale Geschwindigkeit: u_loc = gamma * (e0 + v1*e1 + v2*e2 + v3*e3)
    V4D u_loc = V4D(gamma,gamma*vel.x(1),gamma*vel.x(2),gamma*vel.x(3));

    //double pos2[5] = {pos[0],pos[1],pos[2],pos[3]};
    //mSpacetime->calculateMetric(pos2);
    mSpacetime->localToCoord(pos,u_loc,coordVel,lfType);
  }
  else
  {
    assert(mSolver!=NULL);
    mSolver->startCondition(&pos,coordVel);
  }
  //cerr << "coordvel:\n";
  //coordVel.print(cerr);


  m4d::Metric* metric = mSpacetime;
  GvsLocalTetrad* lt = new GvsLocalTetrad(metric,pos,coordVel);
  lt->adjustTetrad();
  localTetrad.push_back(lt);
  numPositions = localTetrad.size();

  setStartTetrad();
}

//----------------------------------------------------------------------------
//         setStartTetrad
//----------------------------------------------------------------------------
void
GvsStMotionWorldline :: setStartTetrad ( )
{
  //  cerr << "GvsStMotionWorldline :: setStartTetrad()...\n";
  assert (localTetrad[0]!=NULL);
  V4D e0 = localTetrad[0]->getE(0);
  V4D e1 = localTetrad[0]->getE(1);
  V4D e2 = localTetrad[0]->getE(2);
  V4D e3 = localTetrad[0]->getE(3);
  setStartTetrad(e0,e1,e2,e3,true);
}

void
GvsStMotionWorldline :: setStartTetrad ( const V4D &e0, const V4D &e1, const V4D &e2, const V4D &e3,
                                        bool inCoords, GvsLFType lfType  )
{
  //  cerr << "GvsStMotionWorldline :: setStartTetrad( coords )..." << endl;
  //  if (inCoords)  cerr << "yes\n"; else cerr << "no\n";

  // Es muss bereits eine Tetrade mit Startort angegeben sein.
  assert (localTetrad[0]!=NULL);

  V4D e0k = e0;
  V4D e1k = e1;
  V4D e2k = e2;
  V4D e3k = e3;

  //  e0 wird zunaechst abhaengig von der Startbewegungsrichtung gesetzt
  //  und liegt in Koordinatendarstellung vor
  e0k = localTetrad[0]->getVelocity();

  // Sind die Tetradenvektoren e1-e3 bzgl einer natuerlichen Tetrade gegeben, so muessen
  // sie noch auf Koordinatendarstellung transformiert werden.
  m4d::Metric* metric = localTetrad[0]->getMetric();
  P5D pos = localTetrad[0]->getPosition();

  if (!inCoords)
  {
    // metric->localToCoord(pos,e0,e0k,lfType);
    metric->localToCoord(pos,e1,e1k,lfType);
    metric->localToCoord(pos,e2,e2k,lfType);
    metric->localToCoord(pos,e3,e3k,lfType);
  }


  // anschliessend muss noch orthonormiert werden
  GvsGramSchmidt* gs = new GvsGramSchmidt(metric,pos);
  gs->calculateTetrad(e0k,e1k,e2k,e3k);
  delete gs;

  localTetrad[0]->setTetrad(e0k,e1k,e2k,e3k,true);
  localTetrad[0]->setInCoords(true);
}

//----------------------------------------------------------------------------
//         setLambdaStep
//----------------------------------------------------------------------------
void
GvsStMotionWorldline :: setLambdaStep ( double h )
{
  lambdaStep = h;
}

void
GvsStMotionWorldline :: setMaxNumPoints ( int maxPoints )
{
  maxNumPoints = maxPoints;
}

//----------------------------------------------------------------------------
//         calculateGeodesic
//----------------------------------------------------------------------------
void
GvsStMotionWorldline :: calculateGeneralMotion(GvsGeodesicDir geodDir, double deltaT, double tau0, double deltaTau )
{
  //  cerr << "GvsStMotionWorldline::calculateGeodesicMotion\n";
  assert(mSpacetime != NULL);
  assert(mSolver != NULL);
  assert(maxNumPoints > 2);
  
  GvsLocalTetrad* locT;
  locT->setLocalTime(tau0);

  //  Ist bereits ein Teil der Geodaete berechnet, so wird, je nach Richtung in der
  //  die Geodaete fortgesetzt werden soll, die erste oder letzte Tetrade als
  //  Startwert verwendet.
  switch (geodDir)
  {
    case (forward) :
            // localTetradPtr = localTetrad.end();
            locT = localTetrad[numPositions-1];
      break;

    case (backward) :
            // localTetradPtr = localTetrad.begin();
            locT = localTetrad[0];
      break;
  }

  P5D actualPos = locT->getPosition();
  V4D actualVel = locT->getVelocity();

  V4D* base = new V4D[4];
  for (int i=0;i<4;i++)
{
    base[i] = locT->getE(i);
  }

  double min,max;
  mSolver->getBoundingTime(min,max);

  assert( (geodDir==forward) || (geodDir==backward) );
  if (geodDir==forward) mSolver->setGeodCharacter(timelike,geodDir,gvsDestCoDir);
  else mSolver->setGeodCharacter(timelike,geodDir,gvsDestContraDir);

  actualVel = mSolver->getTimeDir()*actualVel;

  mSolver->setBoundingTime(actualPos.x(0),actualPos.x(0)+geodDir*deltaT);
  //    mSolver->print(cerr);

  int numPoints = 0;
  GvsLocalTetrad* localT = new GvsLocalTetrad[maxNumPoints];

  //TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  GvsGeodBreak breakCond =
    mSolver->calcParTransport(actualPos,actualVel,base,lambdaStep,maxNumPoints,localT,numPoints);

  mSolver->errorMessage(breakCond);

  P5D punkt;
  // Starte von i=1, da i=0 schon als Startwert enthalten ist
  for (int i=1;i<numPoints;i++)
  {
    locT = new GvsLocalTetrad(localT[i]);
    locT->setInCoords(true);

    if (geodDir==forward)
    {
      locT->setLocalTime(tau0 + (double(i))*deltaTau);
      localTetrad.push_back(locT);
    }
    else
    {
      locT->setLocalTime(tau0 - (double(i))*deltaTau);
      localTetrad.push_front(locT);
    }
  }
  numPositions = localTetrad.size();

  delete [] localT;
  delete [] base;
}


//----------------------------------------------------------------------------
//         setMotion
//----------------------------------------------------------------------------
void
GvsStMotionWorldline :: setGeneralMotion (const string& filename, bool haveTau )
{
  double tau = 0.0;
  double x0,x1,x2,x3,x4;
  double v0,v1,v2,v3;
  double e0_0, e0_1, e0_2, e0_3;
  double e1_0, e1_1, e1_2, e1_3;
  double e2_0, e2_1, e2_2, e2_3;
  double e3_0, e3_1, e3_2, e3_3;

  ifstream fin(filename.c_str());
  while (true)
  {
    if ( haveTau )
      fin >> tau;
    fin >> x0 >> x1 >> x2 >> x3 >> x4;
    fin >> v0 >> v1 >> v2 >> v3;
    fin >> e0_0 >> e0_1 >> e0_2 >> e0_3;
    fin >> e1_0 >> e1_1 >> e1_2 >> e1_3;
    fin >> e2_0 >> e2_1 >> e2_2 >> e2_3;
    fin >> e3_0 >> e3_1 >> e3_2 >> e3_3;
    P5D pos(x0,x1,x2,x3);
    V4D vel(v0,v1,v2,v3);
    V4D acc(0.0,0.0,0.0,0.0); //TODO brauch man das? ggf numerisch ausrechnen wie in calcpartransport runge kutta (finite differenz)
    V4D e0(e0_0,e0_1,e0_2,e0_3);
    V4D e1(e1_0,e1_1,e1_2,e1_3);
    V4D e2(e2_0,e2_1,e2_2,e2_3);
    V4D e3(e3_0,e3_1,e3_2,e3_3);

    GvsLocalTetrad* lt = new GvsLocalTetrad;
    lt->setInCoords(true);
    lt->setMetric(mSpacetime);
    lt->setPosition(pos);
    lt->setVelocity(vel);
    lt->setAccel(acc);
    lt->setTetrad(e0,e1,e2,e3);//rechnet auch die inverse bereits aus!
    lt->setLocalTime(tau);

    localTetrad.push_back(lt);
    if ( fin.eof())
      break;
  }
  //?? TODO
  localTetrad.pop_back();

  numPositions = localTetrad.size();
}

