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

#include "Obj/Comp/GvsLocalCompObj.h"
#include "Obj/STMotion/GvsStMotionGeodesic.h"

#include "metric/m4dMetric.h"


GvsLocalCompObj :: GvsLocalCompObj() :
    GvsSceneObj(),
    staticTetrad(NULL) {
    objList  = new GvsObjPtrList();

    mNumObjects = 0;
    mPointsWarning = true;
}


GvsLocalCompObj :: GvsLocalCompObj ( GvsStMotion* motion ) : GvsSceneObj(),
    staticTetrad(NULL) {

    if (stMotion!=NULL)
        delete stMotion;
    stMotion = motion;
    haveMotion = true;

    objList  = new GvsObjPtrList();

    mNumObjects = 0;
    mPointsWarning = true;
}

GvsLocalCompObj :: ~GvsLocalCompObj() {
    if (objList != NULL) {
        objList->clear();
        delete objList;

        objList = NULL;
    }

    staticTetrad = NULL;

    mNumObjects = 0;
    mPointsWarning = true;
}

void
GvsLocalCompObj :: setLocalTetrad ( GvsLocalTetrad* locT )
{
    if (haveMotion)
        stMotion->setLocalTetrad(locT);
    else
    {
        staticTetrad = locT;
        staticTetrad->transformTetrad(true); // transformiere auf Koordinaten
    }
}

GvsLocalTetrad*
GvsLocalCompObj :: getLocalTetrad ( int k ) const
{
    std::cerr << "GvsLocalCompObj :: getLocalTetrad\n";
    if (haveMotion)
    {
        assert ( k < numPositions );
        return stMotion->getLocalTetrad(k);
    }
    else
        return staticTetrad;
}


//----------------------------------------------------------------------------
//         set/get motion
//----------------------------------------------------------------------------
void
GvsLocalCompObj :: setMotion ( GvsStMotion *motion )
{
    if (stMotion!=NULL) delete stMotion;
    stMotion = motion;
    haveMotion = true;

    numPositions = stMotion->getNumPositions();

    if (objList->length() != 0) calcSTBoundBoxComplete();
}

GvsStMotion*
GvsLocalCompObj :: getMotion ( void ) const
{
    return stMotion;
}


void GvsLocalCompObj :: Add ( GvsSceneObj *obj ) {
    assert ( obj != NULL );

    objList->Add(obj);
    compBoundBox += obj->boundingBox();
    mNumObjects = objList->length();
    calcSTBoundBox(0);
}


GvsSceneObj* GvsLocalCompObj::getObj ( int nr ) const {
    assert (nr < mNumObjects);
    return objList->getObj(nr);
}

int GvsLocalCompObj::getNumObjs ( void ) const {
    return mNumObjects;
}

/*void
        calcDist  (m4d::Metric* cMetric, const m4d::vec4 &p1, const m4d::vec4 &p2, double &spaceDist, double &timeDist )
        {
          //  p1.print();
          //  p2.print();
          //assert(p1.x(4)==p2.x(4));


          bool cr[4];
          cr[0] =true;
          cr[1] = cr[2] = cr[3] = false;
          cMetric->calculateMetric(p1);

          spaceDist = 0.0;
          timeDist  = 0.0;
          for (int i=0; i<4; i++)
          {
            for (int j=0; j<4; j++)
            {
              if ( !(cr[i]) && !(cr[j]) ) // CoordRange??
               // spaceDist += cMetric->getMetricCoeff(i,j)*(p2.x(i)-p1.x(i))*(p2.x(j)-p1.x(j));
                spaceDist += cMetric->getMetricCoeff(i,j)*cMetric->coordDiff(i,p1.x(i),p2.x(i))*cMetric->coordDiff(j,p1.x(j),p2.x(j));
              if (  (cr[i]) &&  (cr[j]) )
                //timeDist  -= cMetric->getMetricCoeff(i,j)*(p2.x(i)-p1.x(i))*(p2.x(j)-p1.x(j));
                timeDist  -= cMetric->getMetricCoeff(i,j)*cMetric->coordDiff(i,p1.x(i),p2.x(i))*cMetric->coordDiff(j,p1.x(j),p2.x(j));
            }
          }
          if (spaceDist>=0.0) spaceDist = sqrt(spaceDist);
          //else { std::cerr << "SpaceDist < 0 \n"; }
          if (timeDist >=0.0) timeDist  = sqrt(timeDist );
          //else { std::cerr << "TimeDist < 0 \n"; }
        }
*/


bool GvsLocalCompObj :: testIntersection( GvsRay &ray ) {
    // std::cerr << "GvsLocalCompObj :: testIntersection()\n";
    assert (mNumObjects > 0);

    GvsSceneObj* obj = NULL;

    m4d::Metric* stMetric = NULL;
    GvsLocalTetrad* locT0 = NULL;
    GvsLocalTetrad* locT1 = NULL;
    double          localTime0, localTime1;
    m4d::vec4 pos;
    m4d::vec4 box0,box1;

    double t0,t1;
    double spaceDist0,timeDist0,spaceDist1,timeDist1;

    bool result;
    bool intersecFound = false;

    if (!haveMotion) {
        // locT0    = stMotion->getLocalTetrad(0);
        locT0    = staticTetrad;
        stMetric = locT0->getMetric();       // Metrik am Ort der lokalen Tetrade
        pos      = locT0->getPosition();     // Position der lokalen Tetrade

        // box0  = (stMotion->getSTBoundBox())->uppBounds();
        box0 = (locT0->getSTBoundBox())->uppBounds();

        locT1 = locT0;

        int maxSeg   = ray.getNumPoints()-2;
        int startSeg = int(0);
        int endSeg   = maxSeg;

        for (int seg = startSeg; seg < endSeg; seg++) {
            m4d::vec4 p0 = ray.getPoint(seg);
            m4d::vec4 p1 = ray.getPoint(seg+1);

            stMetric->calcSepDist(pos,p0,spaceDist0,timeDist0);
            stMetric->calcSepDist(pos,p1,spaceDist1,timeDist1);

            if ( ((fabs(spaceDist0)<box0.x(1)) && (fabs(timeDist0)<box0.x(0))) ||
                 ((fabs(spaceDist1)<box0.x(1)) && (fabs(timeDist1)<box0.x(0))) )
            {
                // std::cerr << "Hit STBoundBox\n";

                m4d::vec4 p0loc = locT0->transToLocTetrad(p0);
                m4d::vec4 p1loc = locT0->transToLocTetrad(p1);

                for(int i = 0; i < (objList->length()); i++ ) {
                    obj = objList->getObj(i);
                    result = obj->testLocalIntersection(ray,seg,locT0,locT1,p0loc,p1loc,localTime0,localTime1);
                    intersecFound = intersecFound || result;
                }
            }
        }
    }
    else {
        // local object in motion
        assert(stMotion!=NULL);

        int maxSeg   = ray.getNumPoints()-2;
        int startSeg = int(0);
        int endSeg   = maxSeg;

        int numMotionPos = stMotion->getNumPositions();

        for(int seg = startSeg; seg < endSeg; seg++) {
            m4d::vec4 p0 = ray.getPoint(seg);
            m4d::vec4 p1 = ray.getPoint(seg+1);

            t0 = p0.x(0);
            t1 = p1.x(0);

            //fprintf(stderr,"%8.4f %8.4f %8.4f  %8.4f %8.4f.. ",t0,t1,t0*t1,p0.x(1),ray.getTangente(seg).x(0));

            // --------------------
            //  Teste, ob sich der Lichtstrahl und das lokale Objekt zeitlich nahe kommen:
            // ---------------------
            // lokale Tetrade, welche dem Zeitpunkt t0 am naechsten ist
            int num0;
            locT0 = stMotion->getClosestLT(t0,num0);
            if (locT0 == NULL) {
                continue;
            }
            localTime0 = stMotion->getLocalTime(num0);

            stMetric = locT0->getMetric();
            pos      = locT0->getPosition();

            stMetric->calcSepDist( pos,p0,spaceDist0,timeDist0);
            box0  = (locT0->getSTBoundBox())->uppBounds();

            int num1;
            locT1 = stMotion->getClosestLT(t1,num1);
            if (locT1 == NULL) {
                continue;
            }
            localTime1 = stMotion->getLocalTime(num1);

            if ((num0==num1) && (num0!=0))
            {
                if ((num0+1) >= numMotionPos)
                {
                    if ( mPointsWarning )
                    {
                        std::cerr << "Motion has not enough points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
                        mPointsWarning = false;
                        continue;
                    }
                    localTime1 = localTime0;
                }
                else
                {
                    locT1 = stMotion->getLocalTetrad(num0+1);
                    localTime1 = stMotion->getLocalTime(num0+1);
                }
            }

            stMetric = locT1->getMetric();
            pos      = locT1->getPosition();

            stMetric->calcSepDist( pos,p1,spaceDist1,timeDist1);
            box1  = (locT1->getSTBoundBox())->uppBounds();

            //cerr << "box: " << box0.x(1) << " " << box1.x(1) << "  " << box0.x(0) << " " << box1.x(0) << std::endl;
            //cerr << spaceDist0 << " " << spaceDist1 << "  " << timeDist0 << " " << timeDist1 << std::endl;
            //fprintf(stderr,"%8.4f %8.4f %8.4f %8.4f\n",spaceDist0,spaceDist1,timeDist0,timeDist1);

            //printf("--------------  localTime0 localTime1 %f %f\n",localTime0,localTime1);

            // Is one of the points inside the space-time bubble ?
            if ( ((fabs(spaceDist0)<box0.x(1)) && (fabs(timeDist0)<box0.x(0))) ||
                 ((fabs(spaceDist1)<box1.x(1)) && (fabs(timeDist1)<box1.x(0))) )
            {
                // transform points into local tetrad
                m4d::vec4 p0loc = locT0->transToLocTetrad(p0);
                m4d::vec4 p1loc = locT1->transToLocTetrad(p1);

                for(int i = 0; i < (objList->length()); i++ ) {
                    obj = objList->getObj(i);
                    result = obj->testLocalIntersection(ray,seg,locT0,locT1,p0loc,p1loc,localTime0,localTime1);
                    intersecFound = intersecFound || result;
                }
            }
        }
    }
    return intersecFound;
}


GvsBoundBox GvsLocalCompObj::boundingBox( ) const {
    return compBoundBox;
}


void GvsLocalCompObj::calcSTBoundBox ( int k ) {
    m4d::vec3 low = compBoundBox.lowBounds();
    m4d::vec3 upp = compBoundBox.uppBounds();

    double m1 = GVS_MAX(fabs(low.x(0)),fabs(upp.x(0)));
    double m2 = GVS_MAX(fabs(low.x(1)),fabs(upp.x(1)));
    double m3 = GVS_MAX(fabs(low.x(2)),fabs(upp.x(2)));

    GvsLocalTetrad* locT;
    if (haveMotion)
        locT = stMotion->getLocalTetrad(k);
    else
        locT = staticTetrad;

    m4d::vec4 maxBound;
    if (!haveMotion) {
        // an object at rest has a 'time'-cylinder as space-time bubble.
        maxBound = m4d::vec4(DBL_MAX,m1,m2,m3);
    }
    else {
        // ??????????????????????????
        maxBound = m4d::vec4(2.0,m1,m2,m3);
    }
    // Transformiere den Ortsvektor ins Koordinatensystem
    m4d::vec4 maxBoundCoords = locT->transToCoords(maxBound);

    // Bestimme daraus den maximalen raumartigen und zeitartigen Abstand der BoundingBox
    m4d::Metric* metric = locT->getMetric();
    double spaceDist,timeDist;

    //  std::cerr << std::endl;
    //  (locT->getPosition()).print(cerr);
    //  maxBoundCoords.print(cerr);
    metric->calcSepDist(locT->getPosition(),maxBoundCoords,spaceDist,timeDist);
    // std::cerr << spaceDist << " " << timeDist << std::endl;
    // Setze hier spaceDist und timeDist als stBoundingBox
    GvsBoundBox4D* stBoundBox = new GvsBoundBox4D(m4d::vec4(-timeDist,-spaceDist,0,0),m4d::vec4(timeDist,spaceDist,0,0));
    if (haveMotion)
        stMotion->setSTBoundBox(stBoundBox,k);
    else
        staticTetrad->setSTBoundBox(stBoundBox);
    delete stBoundBox;
}


void GvsLocalCompObj::calcSTBoundBoxComplete() {
    // Berechne zu jedem Ort der Bewegung eine stBoundBox und uebergib diese an die zugehoerige
    // lokale Tetrade in der Bewegungsklasse
    //  std::cerr << "GvsLocalCompObj::calcStBoundBoxComplete()...";

    if (haveMotion) {
        int num = stMotion->getNumPositions();
        //  std::cerr << "GvsLocalCompObj :: calcSTBoundBoxComplete: " << num << std::endl;
        for ( int i=0; i < num; i++) {
            calcSTBoundBox(i);
        }
    }
    else {
        calcSTBoundBox();
    }
}


void GvsLocalCompObj::calcSTBoundBoxPartial( int fromPos, int toPos ) {
    if (haveMotion)
    {
        int num = stMotion->getNumPositions();
        assert ( (fromPos < num) && (toPos < num) );

        if (fromPos==toPos) calcSTBoundBox(fromPos);
        else
        {
            int i,von,bis;
            if ( fromPos < toPos )
            {
                von = fromPos;
                bis = toPos;
            }
            else
            {
                von = toPos;
                bis = fromPos;
            }

            i = von;
            do
            {
                calcSTBoundBox(i++);
            }
            while (i<bis);
        }
    }
    else
        calcSTBoundBox();
}


void GvsLocalCompObj :: Print( FILE* fptr ) {
    fprintf(fptr,"LocalCompObject {\n");

    //GvsSceneObj* obj;
    if (!haveMotion) {
        staticTetrad->Print(fptr);
    } else {
        if (stMotion->getMotionType() == gvsMotionGeodesic ) {

        }
        fprintf(fptr,"first point: "); stMotion->getFirstPos().printS(fptr);
        fprintf(fptr,"last point:  "); stMotion->getLastPos().printS(fptr);
    }
    fprintf(fptr,"#objects: %d\n",mNumObjects);
    fprintf(fptr,"}\n");
    /*

  os << "-----------------------------------------\n";
  os << "            LocalCompObject   :" << std::endl;
  os << "-----------------------------------------\n";
  if (!haveMotion)
  {
    os << "Lokale Tetrade (ruhend):" << std::endl;
    staticTetrad->print(os);
  }
  else
  {
    if (stMotion->getMotionType() == gvsMotionGeodesic )
    {
      os << "Geodaetische Bewegung:" << std::endl;
      GvsStMotionGeodesic* gmo = dynamic_cast<GvsStMotionGeodesic*>(stMotion);
      (gmo->getSolver())->print(os);
    }
    os << "Bewegung laeuft vom Punkt: ";
    (stMotion->getFirstPos()).print(os);
    os << "            bis zum Punkt: ";
    (stMotion->getLastPos()).print(os);
  }
  os << "\n\tAnzahl Objekte: " << mNumObjects << std::endl;
  for (int i=0; i<mNumObjects; i++)
  {
    obj = getObj(i);
    obj->print(os);
  }
  compBoundBox.print(os);
  os << "-----------------------------------------\n";
  */
}

