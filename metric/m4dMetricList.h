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
   m4dMetricList.h
 
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

#ifndef  M4D_METRIC_LIST_H
#define  M4D_METRIC_LIST_H

#include "m4dMetric.h"

#include "m4dMetricAlcubierre.h"
#include "m4dMetricAlcubierreSimple.h"
#include "m4dMetricBarriolaVilenkin.h"
#include "m4dMetricBertottiKasner.h"
#include "m4dMetricBesselGravWaveCart.h"
#include "m4dMetricChazyCurzonRot.h"
#include "m4dMetricCosmicStringSchwarzschild.h"
#include "m4dMetricCurzon.h"
#include "m4dMetricDeSitterUniv.h"
#include "m4dMetricDeSitterUnivConf.h"
#include "m4dMetricEddFinkIn.h"
#include "m4dMetricEinsteinRosenWaveWWB.h"
#include "m4dMetricErezRosenVar.h"
#include "m4dMetricErnst.h"
#include "m4dMetricExtremeReissnerNordstromDihole.h"
#include "m4dMetricFriedmanNonEmptyNull.h"
#include "m4dMetricGoedel.h"
#include "m4dMetricGoedelCart.h"
#include "m4dMetricGoedelScaled.h"
#include "m4dMetricGoedelScaledCart.h"
#include "m4dMetricHalilsoyWave.h"
#include "m4dMetricJaNeWi.h"
#include "m4dMetricKasner.h"
#include "m4dMetricKastorTraschen.h"
#include "m4dMetricKerrBL.h"
#include "m4dMetricKottler.h"
#include "m4dMetricMinkowski.h"
#include "m4dMetricMinkowskiConformal.h"
#include "m4dMetricMinkRotLattice.h"
#include "m4dMetricMorrisThorne.h"
#include "m4dMetricPainleveGullstrand.h"
#include "m4dMetricPlaneGravWave.h"
#include "m4dMetricReissnerNordstrom.h"
#include "m4dMetricRotDihole.h"
#include "m4dMetricSchwarzschild.h"
#include "m4dMetricSchwarzschildCart.h"
#include "m4dMetricSchwarzschildIsotropic.h"
#include "m4dMetricSchwarzschildTortoise.h"
#include "m4dMetricStraightSpinningString.h"
#include "m4dMetricSultanaDyer.h"
#include "m4dMetricTaubNUT.h"
#include "m4dMetricTeoWHl.h"
#include "m4dMetricPTD_AI.h"
#include "m4dMetricPTD_AII.h"
#include "m4dMetricPTD_AIII.h"
#include "m4dMetricPTD_BI.h"
#include "m4dMetricPTD_BII.h"
#include "m4dMetricPTD_BIII.h"
#include "m4dMetricPTD_C.h"
#include "m4dMetricPravda_C.h"
#include "m4dMetricPravda_C_Can.h"

namespace m4d
{

const int  NUM_METRICS = 52;

/* --------------------------------------------------------
 *   List of all metrics currently implemented
 *
 *   When editing this list please take care 
 *   of the ordering.
 *
 *   The names here must be equal to the names
 *   given in the constructor of the metric: mMetricName
 * -------------------------------------------------------- */
static const char stl_metric_names[NUM_METRICS][60] =
  {  "unknown",
     "Minkowski",
     "MinkowskiConformal",
     "MinkowskiRotLattice",
     "Schwarzschild",
     "SchwarzschildCart",
     "SchwarzschildIsotropic",
     "SchwarzschildTortoise",
     "EddFinkIn",
     "PainleveGullstrand",
     "AlcubierreWarp",
     "AlcubierreWarpSimple",
     "BarriolaVilenkin",
     "BertottiKasner",
     "BesselGravWaveCart",
     "ChazyCurzonRot",
     "CosmicStringSchwarzschild",
     "Curzon",
     "EinsteinRosenWaveWWB",
     "ErezRosenVar",
     "Ernst",
     "ExtremeReissnerNordstromDihole",
     "FriedmanNonEmptyNull",
     "Goedel",
     "GoedelCart",
     "GoedelScaled",
     "GoedelScaledCart",
     "HalilsoyWave",
     "JanisNewmanWinicour",
     "Kasner",
     "KastorTraschen",
     "KerrBL",
     "Kottler",
     "MorrisThorne",
     "Petrov_Type_D_AI_ES",
     "Petrov_Type_D_AII_ES",
     "Petrov_Type_D_AIII_ES",
     "Petrov_Type_D_BI_ES",
     "Petrov_Type_D_BII_ES",
     "Petrov_Type_D_BIII_ES",
     "Petrov_Type_D_C_ES",
     "PlaneGravWave",
     "Pravda_C-Metric",
     "Pravda_C-Metric_Canonical_Coords",
     "ReissnerNordstrom",
     "RotDihole",
     "DeSitterUniv",
     "DeSitterUnivConformal",
     "StraightSpinningString",
     "SultanaDyerBlackhole",
     "TaubNUT",
     "TeoWHl"
 };

enum  enum_metric
{
    enum_metric_unknown = 0,
    enum_metric_minkowski = 1,
    enum_metric_minkowski_conf,
    enum_metric_minkowski_rotlattice,
    enum_metric_schwarzschild,
    enum_metric_schwarzschild_cart,
    enum_metric_schwarzschild_isotropic,
    enum_metric_schwarzschild_tortoise,
    enum_metric_eddfinkin,
    enum_metric_painleve,
    enum_metric_alcubierre,
    enum_metric_alcubierre_simple,
    enum_metric_barriola,
    enum_metric_bertottikasner,
    enum_metric_bessel_grav_wave_cart,
    enum_metric_chazy_curzon_rot,
    enum_metric_cosmic_string_schwarzschild,
    enum_metric_curzon,
    enum_metric_einstein_rosen_wave_wwb,
    enum_metric_erezrosenvar,
    enum_metric_ernst,
    enum_metric_extreme_reissner_dihole,
    enum_metric_friedman_nonempty_null,
    enum_metric_geodel,
    enum_metric_geodel_cart,
    enum_metric_goedelscaled,
    enum_metric_goedelscaled_cart,
    enum_metric_halilsoy_wave,
    enum_metric_janewi,
    enum_metric_kasner,
    enum_metric_kastortraschen,
    enum_metric_kerrbl,
    enum_metric_kottler,
    enum_metric_morristhorne,
    enum_metric_Petrov_TD_AI,
    enum_metric_Petrov_TD_AII,
    enum_metric_Petrov_TD_AIII,
    enum_metric_Petrov_TD_BI,
    enum_metric_Petrov_TD_BII,
    enum_metric_Petrov_TD_BIII,
    enum_metric_Petrov_TD_C,
    enum_metric_plane_grav_wave,
    enum_metric_Pravda_C_Metric,
    enum_metric_Pravda_C_Can,
    enum_metric_reissner,
    enum_metric_rotdihole,
    enum_metric_desitter_univ,
    enum_metric_desitter_univ_conf,
    enum_metric_spinning_string,
    enum_metric_sultana_dyer,
    enum_metric_taub_nut,
    enum_metric_teowhl
};

} // end namespace m4d

#endif

