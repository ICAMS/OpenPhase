/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Marvin Tegeler
 *
 */

#ifndef MACROS_H
#define MACROS_H

#include <cfenv>

#undef interface

#define OP_LOOP_INTERIOR 0

#define OP_LOOP_WHOLE 1

#ifdef _OPENMP
    #include <omp.h>
#    define myclock_t double
#    define mygettime() omp_get_wtime()
#    define OP_CLOCKS_PER_SEC 1.0
#    define OMP_COLLAPSE_LOOPS 3
#    define OMP_CHUNKSIZE 128
#    define OMP_SCHEDULING_TYPE dynamic
#else
#    define myclock_t clock_t
#    define mygettime() clock()
#    define OP_CLOCKS_PER_SEC CLOCKS_PER_SEC
#endif

#define STORAGE_LOOP_BEGIN(i,j,k,T__,op_loop_bcells__) \
    {\
    if ((long int) op_loop_bcells__ > (T__).Bcells() )\
    {\
        std::cerr << "STORAGE_LOOP: BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_loop_bcells__ << std::endl;\
        throw std::invalid_argument("Bcells");\
    }\
    const long int op_loop_bcells_X__ = std::min((T__).BcellsX(), (long int) op_loop_bcells__); \
    const long int op_loop_bcells_Y__ = std::min((T__).BcellsY(), (long int) op_loop_bcells__); \
    const long int op_loop_bcells_Z__ = std::min((T__).BcellsZ(), (long int) op_loop_bcells__); \
    const long int op_loop_lower_X__ = std::min(-(op_loop_bcells_X__),(long int)0); \
    const long int op_loop_lower_Y__ = std::min(-(op_loop_bcells_Y__),(long int)0); \
    const long int op_loop_lower_Z__ = std::min(-(op_loop_bcells_Z__),(long int)0); \
    const long int op_loop_upper_X__ = std::max((long int)((T__).sizeX() + (op_loop_bcells_X__)),(long int)((T__).sizeX())); \
    const long int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() + (op_loop_bcells_Y__)),(long int)((T__).sizeY())); \
    const long int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() + (op_loop_bcells_Z__)),(long int)((T__).sizeZ())); \
    for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
    for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
    for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
    {

#define STORAGE_LOOP_END \
    } \
}

#define STRINGIFY(...) #__VA_ARGS__

#define OMP_PARALLEL_FOR_LOOP_BEGIN(i,begin__,end__,...) \
{\
    _Pragma(STRINGIFY(omp parallel for schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE) __VA_ARGS__) ) \
    for (long int i = begin__; i < end__; ++i) \
    {

#define OMP_PARALLEL_FOR_LOOP_END \
    } \
}

#define OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,T__,op_loop_bcells__,...) \
{\
    if ((long int) op_loop_bcells__ > (T__).Bcells() )\
    {\
        std::cerr << "OMP_PARALLEL_STORAGE_LOOP: BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_loop_bcells__ << std::endl;\
        throw std::invalid_argument("Bcells");\
    }\
    const long int op_loop_bcells_X__ = std::min((T__).BcellsX(), (long int) op_loop_bcells__); \
    const long int op_loop_bcells_Y__ = std::min((T__).BcellsY(), (long int) op_loop_bcells__); \
    const long int op_loop_bcells_Z__ = std::min((T__).BcellsZ(), (long int) op_loop_bcells__); \
    const long int op_loop_lower_X__ = std::min(-(op_loop_bcells_X__),(long int)0); \
    const long int op_loop_lower_Y__ = std::min(-(op_loop_bcells_Y__),(long int)0); \
    const long int op_loop_lower_Z__ = std::min(-(op_loop_bcells_Z__),(long int)0); \
    const long int op_loop_upper_X__ = std::max((long int)((T__).sizeX() + (op_loop_bcells_X__)),(long int)((T__).sizeX())); \
    const long int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() + (op_loop_bcells_Y__)),(long int)((T__).sizeY())); \
    const long int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() + (op_loop_bcells_Z__)),(long int)((T__).sizeZ())); \
    _Pragma(STRINGIFY(omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE) __VA_ARGS__) ) \
    for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
    for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
    for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
    {

#define OMP_PARALLEL_STORAGE_LOOP_END \
    }\
}
#endif //MACROS_H
