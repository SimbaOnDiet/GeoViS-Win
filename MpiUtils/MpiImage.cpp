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

#include "MpiImage.h"
#include "Img/GvsChannelImg2D.h"
#include "Img/GvsPicIOEnvelope.h"
#include "Img/GvsPictureIO.h"

#include <iostream>
#include <sstream>


GvsMpiImage :: GvsMpiImage()
    : mNumChannels(3),
      mImageWidth(0),
      mImageHeight(0),
      mBuffer(NULL),
      mData(NULL),
      mWithData(false),
      mNumTasksLeft(0),
      mActive(false) {
    mOutfilename = "";
}


GvsMpiImage :: ~GvsMpiImage () {
    if (mBuffer != NULL) {
        delete [] mBuffer;
    }
    if (mData != NULL) {
        delete [] mData;
    }
}

void GvsMpiImage::setWithData( bool withData ) {
    mWithData = withData;
}


void GvsMpiImage::setOutfilename ( const std::string filename, const int num ) {
    std::ostringstream buf;

    std::string fName;
    std::string fExt;
    GvsPictureIO::get_extension(filename, fName, fExt);

    buf << fName << "_" << num << "." << fExt;
    mOutfilename = buf.str();
}


std::string GvsMpiImage :: getOutfilename ( void ) const {
    return mOutfilename;
}


void GvsMpiImage :: setImageResolution ( int x, int y ) {
    assert ((x > 0) && (y > 0));
    mImageWidth  = x;
    mImageHeight = y;
}


void GvsMpiImage :: setNumChannels ( int num ) {
    assert ((num > 0) && (num < 5));
    mNumChannels = num;
}


int GvsMpiImage :: getNumChannels ( void ) const {
    return mNumChannels;
}


void GvsMpiImage :: insertRegion ( int x1, int y1, int x2, int y2, uchar* p, gvsData* data ) {
    if (mActive == false) {
        activate();
    }

    assert (mBuffer != NULL);
    assert (!( (x1 >= mImageWidth) || (x2 >= mImageWidth) || (x1 < 0) || (x2 < 0) ||
               (y1 >= mImageHeight) || (y2 >= mImageHeight) || (y1 < 0) || (y2 < 0)) );

    gvsData* dptr = data;
    for (int j = y1; j <= y2; j++) {
        for (int i = x1; i <= x2; i++) {
            for (int c = 0; c < mNumChannels; c++) {
                mBuffer[(i+mImageWidth*j)*mNumChannels+c] += *p;
                p++;
            }
            if (data!=NULL && mData!=NULL) {
                memcpy(&mData[j*mImageWidth+i],dptr,sizeof(gvsData)*1);
                dptr++;
            }
        }
    }
    decreaseNumTasksLeft();
}


void GvsMpiImage :: setNumTasks ( int num ) {
    mNumTasksLeft = num;
}


int GvsMpiImage :: getNumTasksLeft( void ) const {
    return mNumTasksLeft;
}


void GvsMpiImage :: writePicture( GvsCamFilter filter, double gamma ) {
    std::cerr << "GvsMpiImage :: writePicture\n";
    assert (mOutfilename != "");

    GvsChannelImg2D image(mImageWidth, mImageHeight, mNumChannels, mWithData);
    correctGamma(gamma);
    image.setBlock ( 0, 0, mImageWidth, mImageHeight, mBuffer );
    if (mWithData) {
        image.setDataBlock( 0, 0, mImageWidth, mImageHeight, mData );
        image.writeIntersecData(mOutfilename.c_str(),filter);
    }
    GvsPicIOEnvelope().writeChannelImg(image, (char*)mOutfilename.c_str());
}


bool GvsMpiImage :: writeImageFileIfPossible( GvsCamFilter filter, double gamma ) {
    if (mNumTasksLeft <= 0) {
        writePicture(filter, gamma);
        stop();  // leert den Speicher
        return true;
    }
    return false;
}


void GvsMpiImage :: Print ( FILE* fptr ) const {
    fprintf(fptr,"GvsMpiImage: {\n");
    fprintf(fptr,"\tout:%s  tasks left: %d\n",mOutfilename.c_str(),mNumTasksLeft);
    fprintf(fptr,"}\n");
}


void GvsMpiImage :: correctGamma ( double gamma ) {
    if ( gamma == 0.0 || gamma == 1.0 )
        return;

    double invGamma = 1.0 / gamma;

    int bufferSize = mImageWidth * mImageHeight * mNumChannels;
    double buf;
    for (int i = 0; i < bufferSize; i++)
    {
        buf = pow((double)mBuffer[i], invGamma);
        mBuffer[i] = (uchar)buf;
    }
}


// If there is no buffer yet, do allocate it.
void GvsMpiImage::activate() {
    if (mBuffer == NULL) {
        int bufferSize =  mImageWidth * mImageHeight * mNumChannels;
        mBuffer = new uchar[bufferSize];

        for (int i = 0; i < bufferSize; i++) {
            mBuffer[i] = 0;
        }
    }

    if (mWithData && mData==NULL) {
        int bufferSize = mImageWidth * mImageHeight;
        mData = new gvsData[bufferSize];
    }

    mActive = true;
}


void GvsMpiImage :: decreaseNumTasksLeft ( void ) {
    mNumTasksLeft--;
}


void GvsMpiImage :: stop ( void ) {
    if (mBuffer != NULL) {
        delete [] mBuffer;
    }
    mActive = false;
    mBuffer = NULL;
}
