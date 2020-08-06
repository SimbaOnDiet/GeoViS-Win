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

#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "Img/GvsChannelImg2D.h"
#include "Img/GvsTiffIO.h"

#if HAVE_LIBTIFF==1
#include <tiffio.h>

GvsTiffIO::GvsTiffIO() {
}


bool GvsTiffIO::readChannelImg( GvsChannelImg2D& chanImg, const char * filename) {
    TIFF *tiffImg;
    unsigned int x, y;

    // open the image file
    tiffImg = TIFFOpen( filename,"r");
    if (tiffImg == NULL ) {
        std::cerr << "Error in GvsTiffIO::read(): can't open input file `"
                  << filename << "'." << std::endl;
        return false;
    }

    int numChannels = 3;

    uint32 xsize1, ysize1;
    TIFFGetField(tiffImg, TIFFTAG_IMAGEWIDTH, &xsize1);
    TIFFGetField(tiffImg, TIFFTAG_IMAGELENGTH, &ysize1);
    unsigned int xsize = xsize1;
    unsigned int ysize = ysize1;

    chanImg.setSize(xsize, ysize, numChannels );

    uint32 *buf = new uint32[xsize*ysize];
    assert( buf != NULL );

    int stopOnError = 0;
    if (TIFFReadRGBAImage(tiffImg, xsize, ysize, buf, stopOnError) == 0) {
        std::cerr << "Error in GvsTiffIO::read(): can't convert file `"
                  << filename << "'." << std::endl;
        return false;
    }

    uchar *pixel = new uchar[ numChannels ];
    assert( pixel != NULL );

    //FILE* fptr = fopen("test.ppm","w");
    //fprintf(fptr,"P6\n%d %d\n%d\n",xsize,ysize,255);

    // Copy buffer from Tiff structure to image
    for ( y = 0; y < ysize; y++ ) {
        for ( x = 0; x < xsize; x++ ) {
            pixel[0] = TIFFGetR(buf[x+xsize*y]);  // Red
            pixel[1] = TIFFGetG(buf[x+xsize*y]);  // Green
            pixel[2] = TIFFGetB(buf[x+xsize*y]);  // Blue
            chanImg.setChannels( x, ysize-1-y, pixel );
            //fwrite(pixel,sizeof(uchar),3,fptr);
        }
    }
    //fclose(fptr);

    delete [] pixel;
    delete [] buf;

    chanImg.setFileName( filename );

    TIFFClose( tiffImg );
    return true;
}

bool GvsTiffIO::writeChannelImg( GvsChannelImg2D& chanImg, const char *filename ) {
    if ( chanImg.numChannels() != 3 ) {
        std::cerr << "Error: GvsTiffIO::write() supports only 3-channel-images." << std::endl;
        return false;
    }

    TIFF *out = TIFFOpen(filename, "w");
    if (out == NULL) {
        std::cerr << "Error: GvsTiffIO::write() cannot open file: " << filename << std::endl;
        return false;
    }

    unsigned int nl = static_cast<unsigned int>(chanImg.height());
    unsigned int nc = static_cast<unsigned int>(chanImg.width());

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (uint32) nc);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, (uint32) nl);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
    //TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    // LZW scanline encoding is no longer implemented due to Unisys patent enforcement

    TIFFSetField(out, TIFFTAG_COMPRESSION,1);

    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(out, TIFFTAG_MINSAMPLEVALUE, (uint16) 0);
    TIFFSetField(out, TIFFTAG_MAXSAMPLEVALUE, (uint16) 255);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, (uint32) -1));

    /*
    uchar* ip = chanImg.getImagePtr();
    int pos;
    FILE* fptr = fopen("test.ppm","w");
    fprintf(fptr,"P6\n%d %d\n%d\n",chanImg.width(),chanImg.height(),255);
    for(int i=0; i<chanImg.width()*chanImg.height(); i++) {
        pos = 0*chanImg.width()*chanImg.height() + i;
        fwrite((void*)&ip[pos],sizeof(uchar),1,fptr);
        pos = 1*chanImg.width()*chanImg.height() + i;
        fwrite((void*)&ip[pos],sizeof(uchar),1,fptr);
        pos = 2*chanImg.width()*chanImg.height() + i;
        fwrite((void*)&ip[pos],sizeof(uchar),1,fptr);
    }
    fclose(fptr);
*/
    unsigned int x,y;
    tdata_t buf = _TIFFmalloc(TIFFScanlineSize(out));

    for (y = 0; y < nl; y++) {
        uint8* pp = (uint8*) buf;
        for (x = 0; x < nc; x++) {
            chanImg.sampleChannels( x, y, pp);
            pp += 3;
        }

        if (TIFFWriteScanline(out, buf, y, 0) < 0) {
            _TIFFfree(buf);
            return true;
        }
    }

    _TIFFfree(buf);
    TIFFClose(out);
    return true;
}

#endif
