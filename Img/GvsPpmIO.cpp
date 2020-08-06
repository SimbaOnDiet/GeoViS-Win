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

#include <cstdlib>
#include <iostream>
#include "Img/GvsPpmIO.h"

GvsPpmIO :: GvsPpmIO() {
}

GvsPpmIO :: ~GvsPpmIO() {
}

bool GvsPpmIO :: readChannelImg( GvsChannelImg2D& chanImg, const char * filename) {
//    // cerr << "GvsPpmIO::readChannelImg..." << endl;

    FILE *fptr;
#ifdef _WIN32
    fopen_s(&fptr,filename,"r");
#else
    fptr = fopen(filename,"r");
#endif
    if (fptr==NULL) {
        fprintf(stderr,"Cannot read image %s.\n",filename);
        return false;
    }

    int width,height,c;
    char s[10];
#ifdef _WIN32
    fscanf_s(fptr,"%s\n%d %d\n%d\n",s,&width,&height,&c);
#else
    UNUSED_ATTR int r;
    r = fscanf(fptr,"%s\n%d %d\n%d\n",s,&width,&height,&c);
#endif

    if (width>0 && height>0) {
        char* image;
        if ((image = (char*)malloc(3*width*height*sizeof(char)))==NULL) {
            fprintf(stderr,"Cannot allocate memory for image.\n");
            fclose(fptr);
            return false;
        }
        size_t read = fread(image,sizeof(char),3*width*height,fptr);
        if (read != static_cast<size_t>(3*width*height)) {
            fprintf(stderr,"Possible error reading data!\n");
        }
        fclose(fptr);

        int numChannels = 3;
        chanImg.setSize(width, height, numChannels );
        uchar* col = new uchar[numChannels];

        char* iptr = image;
        for (int i = 0; i < height; i++) {
           for (int j = 0; j < width; j++) {
               col[0] = (uchar)(*iptr++);
               col[1] = (uchar)(*iptr++);
               col[2] = (uchar)(*iptr++);
               chanImg.setChannels(j,i,col);
           }
        }
        delete [] col;
        free(image);
    }
    return true;
}


bool GvsPpmIO :: writeChannelImg( GvsChannelImg2D& chanImg, const char *filename ) {
    // cerr << "GvsPpmIO :: writeChannelImg...\n";

    FILE *fptr;
#ifdef _WIN32
    fopen_s(&fptr,filename,"wb");
#else
    fptr = fopen(filename,"wb");
#endif
    if (fptr==NULL) {
        fprintf(stderr,"Cannot write image %s.\n",filename);
        return false;
    }

    int numChannels = chanImg.numChannels();

    if (numChannels!=3) {
        fprintf(stderr,"NumChannels!=3. Can't write Image as ppm!\n");
        return false;
    }

    int width  = chanImg.width();
    int height = chanImg.height();
    int maxval = 255;
	int c =225 ;


#ifdef _WIN32
   // fscanf_s(fptr,"%s\n%d %d\n%d\n",s,&width,&height,&c); //ORI
   //fscanf_s(fptr,"%s\n%d %d\n%d\n","s",&width,&height,&c);
   fprintf(fptr, "P6\n%d %d\n%d\n", width, height, maxval);
#else
    fprintf(fptr,"P6\n%d %d\n%d\n",width,height,maxval);
#endif

    char* image;
    if ((image = (char*)malloc(3*width*height*sizeof(char)))==NULL) {
        fprintf(stderr,"Cannot allocate memory for image.\n");
        fclose(fptr);
        return false;
    }

    char* iptr = image;
    uchar* col = new uchar[numChannels];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            chanImg.sampleChannels(j,i,col);
            (*iptr++) = col[0];
            (*iptr++) = col[1];
            (*iptr++) = col[2];
        }
    }
    delete [] col;

    size_t written = fwrite(image,sizeof(char),3*width*height,fptr);
    if (written != static_cast<size_t>(3*width*height)) {
        fprintf(stderr,"Possible error writing data!\n");
    }
    fclose(fptr);
    free(image);
    return true;
}
