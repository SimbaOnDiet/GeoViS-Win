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

#include "GvsSolidCSGObj.h"
#include "GvsSolObjSpanList.h"
#include <Ray/GvsRay.h>
#include <Ray/GvsSurfIntersec.h>


GvsSolidCSGObj :: GvsSolidCSGObj( GvsSolidObj *child0,
                                  GvsSolidObj *child1  )
    : GvsSolidObj( NULL )
{
    childSolObj[0] = child0;
    childSolObj[1] = child1;

    if (isValidObjIntersec()) {
        haveObjIntersec = true;
        setChart(child0->getChart());
        setObjType(childSolObj[0]->getObjType());
    }
    else {
        haveObjIntersec = false;
        setChart(0);
        setObjType(local);
    }
}


void GvsSolidCSGObj :: scale( const m4d::vec3& scaleVec ) {
    childSolObj[0]->scale( scaleVec );
    childSolObj[1]->scale( scaleVec );
}

void GvsSolidCSGObj :: translate( const m4d::vec3& transVec ) {
    childSolObj[0]->translate( transVec );
    childSolObj[1]->translate( transVec );
}

void GvsSolidCSGObj :: rotate( const m4d::vec3& rotAxis, double rotAngle ) {
    childSolObj[0]->rotate( rotAxis, rotAngle );
    childSolObj[1]->rotate( rotAxis, rotAngle );
}

void GvsSolidCSGObj :: transform(const m4d::Matrix<double,3,4> &mat ) {
    childSolObj[0]->transform( mat );
    childSolObj[1]->transform( mat );
}


bool GvsSolidCSGObj :: testIntersection( GvsRay& ray ) {
    GvsSolObjSpanList spanList;

    if (haveObjIntersec) {
        if (getRaySpanList( ray, spanList )) {

            GvsSurfIntersec spanStart,spanEnd;
            if ( spanList.getFirstValidSpan( spanStart, spanEnd ) ) {
                if ( ray.isValidSurfIntersec( spanStart.dist() ) ) {
                    spanStart.setSurfIsSelfDescribing( false );
                    //               spanStart.setPoint( ray.origin() + spanStart.dist() * ray.direction());
                    return ray.store( spanStart );
                }
            }
        }
    }
    return false;
}


void GvsSolidCSGObj :: calcNormal( GvsSurfIntersec & ) const {
    std::cerr << "Error in GvsSolidCSGObj::calcNormal(): "
         << "should never be called." << std::endl;
}

void GvsSolidCSGObj :: calcTexUVParam( GvsSurfIntersec & ) const {
    std::cerr << "Error in GvsSolidCSGObj::calcTexUVParam(): "
         << "should never be called." << std::endl;
}

void GvsSolidCSGObj :: calcDerivatives( GvsSurfIntersec & ) const {
    std::cerr << "Error in GvsSolidCSGObj::calcDerivatives(): "
         << "should never be called." << std::endl;
}

void GvsSolidCSGObj :: calcTex3dPoint( GvsSurfIntersec & ) const {
    std::cerr << "Error in GvsSolidCSGObj::calcTex3dPoint(): "
         << "should never be called." << std::endl;
}

void GvsSolidCSGObj :: Print( FILE* ) {
    //std::cerr << "Error in GvsSolidCSGObj::print(ostream& os): "
//         << "should never be called." << std::endl;
}

bool GvsSolidCSGObj :: isValidObjIntersec() {

    //  std::cerr << "objType1 = " << childSolObj[0]->getObjType() << "\t objType2 = "
    //       << childSolObj[1]->getObjType() << std::endl;

    // GvsObjType child1Type = childSolObj[0]->getObjType();
    // GvsObjType child2Type = childSolObj[1]->getObjType();

    //Ueberpfruefe, ob es sich um zwei Objekte gleichen Typs handelt (local, inCoords) und setze das CSG als ein
    //solches Objekt auf, mit der zugehoerigen Karte des ersten der beiden Objekte
    if ( childSolObj[0]->getObjType() == childSolObj[1]->getObjType() )
    {
        if ( childSolObj[0]->getObjType() == inCoords )//Abfrage fuer Koordinatenobjekte
        {
            if ( childSolObj[0]->getChart() != childSolObj[1]->getChart() )
            {
                std::cerr << "GvsSolidCSGObj::isValidObjIntersec(): "
                     << "found coordinate objects are not in the same chart." << std::endl;

                std::cerr << "GvsSolidCSGObj::isValidObjIntersec(): "
                     << "obj1-chart = " << childSolObj[0]->getChart()
                     << "   obj2-chart = " << childSolObj[1]->getChart() << std::endl;


                this->setChart  ( childSolObj[0]->getChart() );
                exit(0);
                return false;
            }
            else
            {
                //            std::cerr << "GvsSolidCSGObj::isValidObjIntersec(): "
                //                 << "found objects were successfully intersected." << std::endl;

                this->setObjType( childSolObj[0]->getObjType() );
                this->setChart  ( childSolObj[0]->getChart()   );   //Muss nur fuer Koordinatenobjekte gemacht
                //werden,da lokale Objekte in der Tetrade
                //leben, und da gibt es nur eine Karte
                //            std::cerr << "CSG-Obj: ObjType: child1 = " << childSolObj[0]->getObjType() << std::endl;
                //            std::cerr << "CSG-Obj: ObjType: child2 = " << childSolObj[1]->getObjType() << std::endl;
                //	    std::cerr << "CSG-Obj: ObjType: Type = " << this->getObjType() << std::endl;
                //	    std::cerr << "CSG-Obj: Chart = " << this->getChart() << std::endl << std::endl;
                return true;
            }
        }

        else //Abfrage fuer lokale Objekte
        {
            //         std::cerr << "GvsSolidCSGObj::isValidObjIntersec(): "
            //              << "found objects were successfully intersected." << std::endl;
            this->setObjType( childSolObj[0]->getObjType() );
            return true;
        }
    }

    else
    {
        std::cerr << "GvsSolidCSGObj::isValidObjIntersec(): "
             << "You cannot intersect local- and inCoord-objects." << std::endl;

        //      std::cerr << "CSG-Obj: ObjType: child1 = " << childSolObj[0]->getObjType() << std::endl;
        //      std::cerr << "CSG-Obj: ObjType: child2 = " << childSolObj[1]->getObjType() << std::endl;
        //      std::cerr << "CSG-Obj: ObjType: Type = " << this->getObjType() << std::endl;
        //      std::cerr << "CSG-Obj: Chart = " << this->getChart() << std::endl << std::endl;
        exit(0);
        return false;
    }
}

GvsCSGType GvsSolidCSGObj :: getCSGType () const {
    return mCSGType;
}

void GvsSolidCSGObj :: setCSGType( GvsCSGType csgType) {
    mCSGType = csgType;
}
