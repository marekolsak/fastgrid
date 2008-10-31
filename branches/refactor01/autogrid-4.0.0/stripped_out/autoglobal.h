/*

 $Id: autoglobal.h,v 1.5 2007/05/03 20:46:06 garrett Exp $

 AutoGrid 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#pragma once
#include <sys/types.h>
#include <cstring>
#include <cstdio>
#include "structs.h"

/******************************************************************************/
/*      Name: autoglobal.h                                                    */
/*  Function: Global variables for Autodock modules.                          */
/* Copyright: (C) TSRI                                                        */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett M. Morris                                               */
/*                                                                            */
/*            e-mail: garrett@scripps.edu				      */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*      Date: 03/18/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Nothing.                                                        */
/*   Globals: programname, AutoGridHelp, AutoDockHelp, command_mode,          */
/*            command_in_fp, command_out_fp, GPF, logFile.                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 03/18/93 GMM     Created.                                                  */
/******************************************************************************/


/*----------------------------------------------------------------------------*/
/* Globals,                                                                   */
/*----------------------------------------------------------------------------*/

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// All global variables were moved.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/*
// AutoDock 3 Linear Free Energy Model Coefficients wrt AD2
AD3.coeff_vdW     = 0.1485L;
AD3.coeff_hbond   = 0.0656L;
AD3.coeff_estat   = 0.1146L;
AD3.coeff_desolv  = 0.1711L;
AD3.coeff_tors    = 0.3113L;

// AutoDock 3 Linear Free Energy Model Standard Errors wrt AD2
AD3.stderr_vdW    = 0.0237L;
AD3.stderr_hbond  = 0.0558L;
AD3.stderr_estat  = 0.0238L;
AD3.stderr_desolv = 0.1035L;
AD3.stderr_tors   = 0.0910L;

// AutoDock 4 Linear Free Energy Model Coefficients wrt AD3
AD4_wrt_3.coeff_vdW    = 1.002L;
AD4_wrt_3.coeff_hbond  = 1.931L;
AD4_wrt_3.coeff_estat  = 1.229L;
AD4_wrt_3.coeff_desolv = 0.122L;
AD4_wrt_3.coeff_tors   = 0.290L;

// AutoDock 4 Linear Free Energy Model Standard Errors wrt AD3
AD4_wrt_3.stderr_vdW    = 0.059L;
AD4_wrt_3.stderr_hbond  = 0.329L;
AD4_wrt_3.stderr_estat  = 0.173L;
AD4_wrt_3.stderr_desolv = 0.026L;
AD4_wrt_3.stderr_tors   = 0.041L;

// AutoDock 4 Linear Free Energy Model Coefficients wrt AD2
AD4.coeff_vdW    = AD4_wrt_3.coeff_vdW    * AD3.coeff_vdW;
AD4.coeff_hbond  = AD4_wrt_3.coeff_hbond  * AD3.coeff_hbond;
AD4.coeff_estat  = AD4_wrt_3.coeff_estat  * AD3.coeff_estat;
AD4.coeff_desolv = AD4_wrt_3.coeff_desolv * 1.0L;
AD4.coeff_tors   = AD4_wrt_3.coeff_tors   * AD3.coeff_tors;

// AD4 FE
//        coeff
// vdW    0.148797
// hbond  0.1266736
// estat  0.1408434
// desolv 0.122
// tors   0.090277
*/
