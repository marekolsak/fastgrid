/*

 $Id: partokens.h,v 1.2 2007/05/03 20:46:06 garrett Exp $

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


/******************************************************************************
 *      Name: partoken.h                                                      *
 *  Function: Define tokens for parsing AutoDock atomic parameter files       *
 * Copyright: (C) Garrett Matthew Morris, TSRI                                *
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: 08/03/2005                                                      * 
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 08/03/05 GMM     Created                                                   *
 ******************************************************************************/

#ifndef PAR_TOKENS
#define PAR_TOKENS

#define	PAR_		-1
#define	PAR_NULL	 0
#define	PAR_VDW 	 1
#define	PAR_HBOND	 2
#define	PAR_ESTAT	 3
#define	PAR_DESOLV	 4
#define	PAR_TORS	 5
#define	PAR_ATOM_PAR	 6
#define	PAR_COMMENT	 7

#endif
/* EOF */
