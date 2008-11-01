/*

 $Id: ranlib.h,v 1.3 2007/05/03 20:46:06 garrett Exp $

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

/* Prototypes for all user accessible RANLIB routines */

#pragma once
#include "typedefs.h"

void advnst(FourByteLong k);
Real genbet(Real aa,Real bb);
Real genchi(Real df);
Real genexp(Real av);
Real genf(Real dfn, Real dfd);
Real gengam(Real a,Real r);
void genmn(Real *parm,Real *x,Real *work);
void genmul(FourByteLong n,Real *p,FourByteLong ncat,FourByteLong *ix);
Real gennch(Real df,Real xnonc);
Real gennf(Real dfn, Real dfd, Real xnonc);
Real gennor(Real av,Real sd);
void genprm(FourByteLong *iarray,int larray);
Real genunf(Real low,Real high);
void getsd(FourByteLong *iseed1,FourByteLong *iseed2);
void gscgn(FourByteLong getset,FourByteLong *g);
FourByteLong ignbin(FourByteLong n,Real pp);
FourByteLong ignnbn(FourByteLong n,Real p);
FourByteLong ignpoi(Real mu);
FourByteLong ignuin(FourByteLong low,FourByteLong high);
void initgn(FourByteLong isdtyp);
FourByteLong mltmod(FourByteLong a,FourByteLong s,FourByteLong m);
void phrtsd(char* phrase,FourByteLong* seed1,FourByteLong* seed2);
void setall(FourByteLong iseed1,FourByteLong iseed2);
void setant(FourByteLong qvalue);
void setgmn(Real *meanv,Real *covm,FourByteLong p,Real *parm);
void setsd(FourByteLong iseed1,FourByteLong iseed2);
Real sgamma(Real a);
Real rcauchy(Real, Real);

extern FourByteLong ignlgi;
extern Real ranf;
extern Real sexpo;
extern Real snorm;
extern Real scauchy1;
