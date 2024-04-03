/*
 * ligand/dSFMT-params19937.h
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef DSFMT_PARAMS19937_H
#define DSFMT_PARAMS19937_H

/* #define DSFMT_N	191 */
/* #define DSFMT_MAXDEGREE	19992 */
#define DSFMT_POS1	117
#define DSFMT_SL1	19
#define DSFMT_MSK1	UINT64_C(0x000ffafffffffb3f)
#define DSFMT_MSK2	UINT64_C(0x000ffdfffc90fffd)
#define DSFMT_MSK32_1	0x000ffaffU
#define DSFMT_MSK32_2	0xfffffb3fU
#define DSFMT_MSK32_3	0x000ffdffU
#define DSFMT_MSK32_4	0xfc90fffdU
#define DSFMT_FIX1	UINT64_C(0x90014964b32f4329)
#define DSFMT_FIX2	UINT64_C(0x3b8d12ac548a7c7a)
#define DSFMT_PCV1	UINT64_C(0x3d84e1ac0dc82880)
#define DSFMT_PCV2	UINT64_C(0x0000000000000001)
#define DSFMT_IDSTR	"dSFMT2-19937:117-19:ffafffffffb3f-ffdfffc90fffd"


/* PARAMETERS FOR ALTIVEC */
#if defined(__APPLE__)	/* For OSX */
    #define ALTI_SL1 	(vector unsigned char)(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
    #define ALTI_SL1_PERM \
	(vector unsigned char)(2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1)
    #define ALTI_SL1_MSK \
	(vector unsigned int)(0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U)
    #define ALTI_MSK	(vector unsigned int)(DSFMT_MSK32_1, \
			DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4)
#else	/* For OTHER OSs(Linux?) */
    #define ALTI_SL1 	{3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}
    #define ALTI_SL1_PERM \
	{2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1}
    #define ALTI_SL1_MSK \
	{0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U}
    #define ALTI_MSK \
	{DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4}
#endif

#endif /* DSFMT_PARAMS19937_H */
