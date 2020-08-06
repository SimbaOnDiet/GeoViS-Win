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
    m4dUtilities.cpp 
 
  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave
 
 
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

#include "m4dUtilities.h"

namespace m4d 
{

// ---------------------------------------------------
/*!  Get system clock in micro-seconds.
 */
int64_t  get_system_clock() 
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (int64_t)tv.tv_sec * 1000000 + tv.tv_usec;
}


// ---------------------------------------------------
/*!  Tokenize file
 *  \param filename : name of file.
 *  \param tokens   : reference to vector of vector of string.
 *  \param useStandardIgnoreTokens : use standard ignore tokens ("#").
 *  \return  true : success.
 */
bool tokenizeFile( const std::string filename, std::vector<std::vector<std::string> > &tokens, bool useStandardIgnoreTokens )
{
  std::ifstream in(filename.c_str());
  
  if (!in) {
    fprintf(stderr,"Cannot open file \"%s\"\n",filename.c_str());
    return false;
  }
    
  do
  {
    std::string line;
    getline(in, line);
  
    bool ignoreLine = false;  
    if (!line.compare(0,1,""))
      ignoreLine = true;
    
    if (useStandardIgnoreTokens)
    {
      if (!line.compare(0,1,"#"))
        ignoreLine = true;
    }
    
    if (!ignoreLine)
    {
      std::string buf;
      std::stringstream ss(line);
    
      std::vector<std::string> line_tokens;
  
      while(ss >> buf)
        line_tokens.push_back(buf);
  
      tokens.push_back(line_tokens);
    }
  } while (!in.eof());
  
  in.close();
  return true;
}

// ---------------------------------------------------
/*!  Tokenize file
 *  \param filename : name of file.
 *  \param ignores  : vector of strings which indicates lines that have to be ignored.
 *  \param tokens   : reference to vector of vector of string.
 *  \return  true : success.
 */
bool tokenizeFile ( const std::string filename, const std::vector<std::string> ignores, std::vector<std::vector<std::string> > &tokens )
{
  std::ifstream in(filename.c_str());
  
  if (!in) {
    fprintf(stderr,"Cannot open file %s\n",filename.c_str());
    return false;
  }
  
  do
  {
    std::string line;
    getline(in, line);
  
    bool ignoreLine = false;
   
    for(unsigned int i=0; i<ignores.size(); i++)
      if (!line.compare(0,ignores[i].length(),ignores[i]))
        ignoreLine |= true;
  
    if (!line.compare(0,1,""))
      ignoreLine = true;
  
    if (!ignoreLine)
    {    
      std::vector<std::string> line_tokens;  
      std::string buf;
      std::stringstream ss(line);
  
      while(ss >> buf)
        line_tokens.push_back(buf);
  
      tokens.push_back(line_tokens);
    }
  } while (!in.eof());

  in.close();
  return true;
}

// ---------------------------------------------------
/*!  get integer from tokens 
 */
bool  getIntFromTokens ( const std::vector<std::string> &tokenRow, std::string name, int &val )
{
   std::string baseString = tokenRow[0];
   if (baseString.compare(name)==0 && tokenRow.size()>1)
   { 
      val = atoi(tokenRow[1].c_str());
      return true;
   }
   return false;
}

// ---------------------------------------------------
/*!  get integer array from tokens 
 */
bool getIntFromTokensV    ( const std::vector<std::string> &tokenRow, std::string name, int num, int* val )
{
   std::string baseString = tokenRow[0];
   if (baseString.compare(name)==0 && (int)tokenRow.size()>(num+1))
   { 
      for(int i=0; i<num; i++)
        val[i] = atoi(tokenRow[i+1].c_str());
      return true;
   }
   return false;  
}

// ---------------------------------------------------
/*!  get double from tokens 
 */
bool getDoubleFromTokens  ( const std::vector<std::string> &tokenRow, std::string name, double &val )
{
   std::string baseString = tokenRow[0];
   if (baseString.compare(name)==0 && tokenRow.size()>1)
   { 
      val = atof(tokenRow[1].c_str());
      return true;
   }
   return false;
}

// ---------------------------------------------------
/*!  get double array from tokens 
 */
bool getDoubleFromTokensV ( const std::vector<std::string> &tokenRow, std::string name, int num, double* val )
{
   std::string baseString = tokenRow[0];
   if (baseString.compare(name)==0 && (int)tokenRow.size()>(num+1))
   { 
      for(int i=0; i<num; i++)
        val[i] = atof(tokenRow[i+1].c_str());
      return true;
   }
   return false;  
}

// ---------------------------------------------------
/*!  get string from tokens 
 */
bool getStringFromTokens  ( const std::vector<std::string> &tokenRow, std::string name, std::string &val )
{
   std::string baseString = tokenRow[0];
   if (baseString.compare(name)==0 && tokenRow.size()>1)
   { 
      val = tokenRow[1];
      return true;
   }
   return false;
}

/*! Write a binary float array.
 *   \param   filename  :  name of the file.
 *   \param      array  :  pointer to float array.
 *   \param          x  :  width of array.
 *   \param          y  :  height of array.
 *   \param          c  :  number of channels.
 *   \return      true  :  success.
 *   \return     false  :  error occured.
 */
bool  writeFloatArray ( std::string filename, const float* array, int x, int y, int c )
{
  std::ofstream out(filename.c_str(), std::ios::binary);
  if (!out.is_open())
  {
     fprintf(stderr,"Can't open file %s for output.\n",filename.c_str());
     return false;
  }

  std::string hdr = "HEAD";
  char buf[5];
#ifdef _WIN32
  sprintf_s(buf,"%s",hdr.c_str());
#else
  sprintf(buf,"%s",hdr.c_str());
#endif
  std::string hdrd = "DATA";
  char bufd[5];
#ifdef _WIN32
  sprintf_s(bufd,"%s",hdrd.c_str());
#else
  sprintf(bufd,"%s",hdrd.c_str());
#endif
  
  int hdrSize = 24;
  out.write((char*)&buf[0], sizeof(char)*4);
  out.write((char*)&hdrSize, sizeof(int));  
  out.write((char*)&x, sizeof(int));
  out.write((char*)&y, sizeof(int));
  out.write((char*)&c, sizeof(int));
  out.write((char*)&bufd[0], sizeof(char)*4);
  out.write((char*)&array[0], sizeof(float)*x*y*c);
  out.close();
  return true;  
}

/*! Read a binary float array.
 *   \param   filename  :  name of the file.
 *   \param          x  :  reference to width of array.
 *   \param          y  :  reference to height of array.
 *   \param          c  :  reference to number of channels.
 *   \return      true  :  success.
 *   \return     false  :  error occured.
 */
float* readFloatArray   ( std::string filename, int &x, int &y, int &c )
{
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in.is_open())
  {
     fprintf(stderr,"Can't open file %s for reading.\n",filename.c_str());
     return false;
  }
  
  char buf[5];
  char bufd[5];
  for(int i=0; i<4; i++)
  {
    buf[i] = in.peek();
    in.seekg(i+1);
  }
  
  float* array = NULL;
  if (strncmp(buf,"HEAD",4)==0)
  {
    int hdrSize;
    in.read((char*)&hdrSize,sizeof(int));   // header size
    in.read((char*)&x,sizeof(int));         // resX
    in.read((char*)&y,sizeof(int));         // resY
    in.read((char*)&c,sizeof(int));         // numChannels
    in.read((char*)&bufd,sizeof(char)*4);
    if (strncmp(bufd,"DATA",4)==0)
    {
      if (array!=NULL)
        delete [] array;
      array = new float[x*y*c];
	  std::streamsize size = x*y*c*sizeof(float);
      in.read((char*)&array[0],size);
    }
    in.close();    
  }
  
  return array;
}

} // end namespace m4d

