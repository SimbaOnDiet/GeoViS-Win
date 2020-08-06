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
/* dynload.c Dynamic Loader for TinyScheme */
/* Original Copyright (c) 1999 Alexander Shendi     */
/* Modifications for NT and dl_* interface, scm_load_ext: D. Souflis */
/* Refurbished by Stephen Gildea */

#define _SCHEME_SOURCE
#include "dynload.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MAXPATHLEN
# define MAXPATHLEN 1024
#endif

static void make_filename(const char *name, char *filename); 
static void make_init_fn(const char *name, char *init_fn); 

#ifdef _WIN32
# include <windows.h>
#else
typedef void *HMODULE;
typedef void (*FARPROC)();
#define SUN_DL
#include <dlfcn.h>
#endif

#ifdef _WIN32

#define PREFIX ""
#define SUFFIX ".dll"

 static void display_w32_error_msg(const char *additional_message) 
 { 
   LPVOID msg_buf; 
 
   FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_ALLOCATE_BUFFER, 
		 NULL, GetLastError(), 0, 
		 (LPTSTR)&msg_buf, 0, NULL); 
   fprintf(stderr, "scheme load-extension: %s: %s", additional_message, msg_buf); 
   LocalFree(msg_buf); 
 } 

static HMODULE dl_attach(const char *module) {
  HMODULE dll = LoadLibrary(module);
  if (!dll) display_w32_error_msg(module); 
  return dll; 
}

static FARPROC dl_proc(HMODULE mo, const char *proc) {
  FARPROC procedure = GetProcAddress(mo,proc); 
  if (!procedure) display_w32_error_msg(proc); 
  return procedure; 
}

static void dl_detach(HMODULE mo) {
 (void)FreeLibrary(mo);
}

#elif defined(SUN_DL)

#include <dlfcn.h>

#define PREFIX "lib"
#define SUFFIX ".so"

static HMODULE dl_attach(const char *module) {
  HMODULE so=dlopen(module,RTLD_LAZY);
  if(!so) {
    fprintf(stderr, "Error loading scheme extension \"%s\": %s\n", module, dlerror()); 
  }
  return so;
}

static FARPROC dl_proc(HMODULE mo, const char *proc) {
  const char *errmsg;
  FARPROC fp=(FARPROC)dlsym(mo,proc);
  if ((errmsg = dlerror()) == 0) {
    return fp;
  }
  fprintf(stderr, "Error initializing scheme module \"%s\": %s\n", proc, errmsg);
 return 0;
}

static void dl_detach(HMODULE mo) {
 (void)dlclose(mo);
}
#endif

pointer scm_load_ext(scheme *sc, pointer args)
{
   pointer first_arg;
   pointer retval;
   char filename[MAXPATHLEN], init_fn[MAXPATHLEN+6];
   char *name;
   HMODULE dll_handle;
   void (*module_init)(scheme *sc);
   
   if ((args != sc->NIL) && is_string((first_arg = pair_car(args)))) {
      name = string_value(first_arg);
      make_filename(name,filename);     
      make_init_fn(name,init_fn);     
      dll_handle = dl_attach(filename);
      if (dll_handle == 0) {
         retval = sc -> F;
      }
      else {
         module_init = (void(*)(scheme *))dl_proc(dll_handle, init_fn);
         if (module_init != 0) {
            (*module_init)(sc);
            retval = sc -> T;
         }
         else {
            retval = sc->F;
         }
      }
   }
   else {
      retval = sc -> F;
   }
   
  return(retval);
}

static void make_filename(const char *name, char *filename) {
 strcpy(filename,name);
 strcat(filename,SUFFIX);
}         

static void make_init_fn(const char *name, char *init_fn) {
 const char *p=strrchr(name,'/');
 if(p==0) {
     p=name;
 } else {
     p++;
 }
 strcpy(init_fn,"init_");
 strcat(init_fn,p);
}






