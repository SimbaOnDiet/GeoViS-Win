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

#include <vector>
#include <fstream>

#include "MpiTaskManager.h"

#include "Img/GvsPicIOEnvelope.h"
#include "Img/GvsChannelImg2D.h"

#include "m4dGlobalDefs.h"

extern std::vector<GvsDevice*>  gpDevice;


GvsMpiTaskManager :: GvsMpiTaskManager( const std::string infile, const std::string outfile ) {
    inFileName  = infile;
    outFileName = outfile;

    parser = new GvsParser();

    mRenderDevice = -1;  // alle Devices Rendern
    mStartDevice = 0;
}


GvsMpiTaskManager::~GvsMpiTaskManager ( ) {
    delete parser;
}


bool GvsMpiTaskManager :: setStartDevice ( int startdev ) {
    mStartDevice = startdev;
    if (mStartDevice<0)
        mStartDevice = 0;
    return true;
}


int GvsMpiTaskManager :: getStartDevNr() const {
    return mStartDevice;
}


void GvsMpiTaskManager :: setRenderDevice ( int renderdev ) {
    mRenderDevice = renderdev;
    if (mRenderDevice!=-1)
        mStartDevice  = renderdev;
}


bool GvsMpiTaskManager :: initialize ( int , int numNodesImage ) {
    //fprintf(stderr,"Parse scm file...\n");
    parser->read_scene(inFileName.c_str());
    mNumDevices = parser->getNumDevices();

    if (mRenderDevice!=-1) {
        if (mRenderDevice >= mNumDevices) {
            fprintf(stderr,"Device %d does not exist!\n",mRenderDevice);
            return false;
        }
        fprintf(stderr,"Render only device %d of %d\n",mRenderDevice,mNumDevices);
        mNumDevices = mRenderDevice+1;
    }

    mImage = new GvsMpiImage[mNumDevices-mStartDevice];
    parser->getDevice(&mDevice,0);
    m4d::ivec2 res(mDevice.camera->GetResolution());

    mImageWidth  = res.x(0);
    mImageHeight = res.x(1);

    mNumNodesImage = numNodesImage;

    // Bildgroesse kleiner als die Anzahl Knoten je Bild
    if (mImageHeight < mNumNodesImage )
        mNumNodesImage = mImageHeight;

    // Schrittweite
    int yStepsize = mImageHeight / mNumNodesImage;
    if (yStepsize < 1)
        yStepsize = 1;

    mNumNodesImage = mImageHeight / yStepsize;
    if (mNumNodesImage*yStepsize < mImageHeight)
        mNumNodesImage++;


    // Anzahl der Tasks ergibt sich aus der Anzahl Teilaufgaben je Bild mal der Anzahl Bilder,
    // wobei die Anzahl der Bilder gleich der Anzahl der Devices ist
    mNumTasks = mNumNodesImage * (mNumDevices-mStartDevice);
    mTasks = new GvsMpiTask[mNumTasks];
    assert( mTasks != NULL );

    GvsCamFilter  camFilter;

    // Subdivide images into tasks
    int actTask = 0;
    for ( int image = 0; image < mNumDevices-mStartDevice; image++)   {
        mImage[image].setNumTasks(mNumNodesImage);
        mImage[image].setImageResolution(mImageWidth,mImageHeight);
        mImage[image].setOutfilename(outFileName,image+mStartDevice);

        camFilter = mDevice.camera->getCamFilter();

        //int numChannels;
        switch (camFilter) {
            default:
            case gvsCamFilterRGB: {
                //numChannels = 3;
                break;
            }
            case gvsCamFilterRGBpdz:
            case gvsCamFilterRGBjac: {
                //numChannels = 3;
                mImage[image].setWithData(true);
                break;
            }
        }

        int startY = 0;
        for ( int imageRegion = 0; imageRegion < mNumNodesImage; imageRegion++ ) {
            mTasks[actTask].x1 = 0;
            mTasks[actTask].y1 = startY;
            mTasks[actTask].x2 = mImageWidth-1;
            startY += yStepsize;
            mTasks[actTask].y2 = startY-1;

            mTasks[actTask].status = TASK_WAITING;

            // assign a device to each task
            parser->getDevice(&(mTasks[actTask].device),image+mStartDevice);

            // Every task stores in specific image.
            mTasks[actTask].imageNr = image;

            actTask++;
        }
        mTasks[actTask-1].y2 = mImageHeight-1;
    }
    return true;
}


void GvsMpiTaskManager :: getDevice  ( GvsDevice *device, unsigned int k ) {
    assert ( k < gpDevice.size() );

    device->metric = gpDevice[k]->metric;

    device->camera    = gpDevice[k]->camera;
    device->projector = gpDevice[k]->projector;

    device->lightSrcMgr = gpDevice[k]->lightSrcMgr;
    device->sceneGraph  = gpDevice[k]->sceneGraph;

    device->mChangeObj = gpDevice[k]->mChangeObj;
    device->listNumChangeObj = gpDevice[k]->listNumChangeObj;
}


int GvsMpiTaskManager::getNumTasks ( ) const {
    return mNumTasks;
}


void GvsMpiTaskManager::activateTask ( int task ) {
    mTasks[task].status = TASK_RUNNING;
}


void GvsMpiTaskManager::getViewPort ( int task, int &x1, int &y1, int &x2, int &y2) const {
    x1 = mTasks[task].x1;
    y1 = mTasks[task].y1;
    x2 = mTasks[task].x2;
    y2 = mTasks[task].y2;
}


int GvsMpiTaskManager::getImageNr ( int task ) const {
    return mTasks[task].imageNr;
}


void GvsMpiTaskManager :: insertRegion ( int task, uchar* p, gvsData* data ) {
    //  cerr << "GvsMpiTaskManager :: insertRegion: " << task << endl;
    int x1 = mTasks[task].x1;
    int y1 = mTasks[task].y1;
    int x2 = mTasks[task].x2;
    int y2 = mTasks[task].y2;

    mImage[mTasks[task].imageNr].insertRegion(x1,y1,x2,y2,p,data);
}


bool GvsMpiTaskManager :: writeImageFileIfPossible ( int task, double gamma ) {
    GvsCamFilter filter = mDevice.camera->getCamFilter();
    return mImage[mTasks[task].imageNr].writeImageFileIfPossible( filter, gamma );
}


int GvsMpiTaskManager :: getNextAvailableTask ( int task ) const {
    while ((mTasks[task].status != TASK_WAITING) && (task < mNumTasks-1))   {
        task++;
    }

    if (task == mNumTasks-1)   {
        if (mTasks[task].status != TASK_WAITING) {
            task++;
        }
    }
    else if (task >= mNumTasks) {
        task = mNumTasks;
    }

    return task;
}


void GvsMpiTaskManager :: createScene ( int task, GvsDevice *device ) {
    device->metric   = (mTasks[task].device).metric;

    device->camera      = (mTasks[task].device).camera;
    device->projector   = (mTasks[task].device).projector;

    device->lightSrcMgr = (mTasks[task].device).lightSrcMgr;
    device->sceneGraph  = (mTasks[task].device).sceneGraph;

    device->mChangeObj  = (mTasks[task].device).mChangeObj;
    device->listNumChangeObj = (mTasks[task].device).listNumChangeObj;

    // make changes
    device->makeChange();
}


void GvsMpiTaskManager :: Print ( FILE* fptr ) const {
    fprintf(fptr,"\nMpiTaskManager:   \n---------------\n");
    fprintf(fptr,"\tno devices:      %d\n",mNumDevices);
    fprintf(fptr,"\tno nodes/image:  %d\n",mNumNodesImage);
    fprintf(fptr,"\tno tasks:        %d\n",mNumTasks);
    fprintf(fptr,"\tscene file:      %s\n",inFileName.c_str());
    fprintf(fptr,"\timage file:      %s\n",outFileName.c_str());
    fprintf(fptr,"\n");
}
