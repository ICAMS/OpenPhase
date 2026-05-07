/*
 *  This file is part of the OpenPhase (R) software library.
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
 *
 *  File created :   2023
 *  Main contributors :   Marvin Tegeler; Oleg Shchyglo;
 *
 */

#ifndef OP_MPI_H
#define OP_MPI_H

#include <sstream>
#include <vector>
#include "fftw3.h"

enum OP_MPI_Datatype
{
	OP_MPI_CHAR,
    OP_MPI_INT,
    OP_MPI_LONG,
    OP_MPI_FLOAT,
    OP_MPI_DOUBLE,
    OP_MPI_DOUBLE_INT,
    OP_MPI_UNSIGNED,
    OP_MPI_UNSIGNED_LONG,
    OP_MPI_UNSIGNED_LONG_LONG,
    OP_MPI_CXX_BOOL
};
enum OP_MPI_Op
{
    OP_MPI_SUM,
    OP_MPI_PROD,
    OP_MPI_MIN,
    OP_MPI_MAX,
    OP_MPI_MAXLOC,
    OP_MPI_LOR,
    OP_MPI_LAND
};
enum OP_MPI_Status
{
    OP_MPI_STATUS_IGNORE,                                                       ///< Dummy replacement for MPI_STATUS_IGNORE
};
enum OP_MPI_Comm
{
    OP_MPI_COMM_WORLD,
    OP_MPI_COMM_SELF
};
enum OP_MPI_Options
{
    OP_MPI_IN_PLACE
};
enum OP_MPI_Threads
{
    OP_MPI_THREAD_SINGLE,
    OP_MPI_THREAD_FUNNELED,
    OP_MPI_THREAD_SERIALIZED,
    OP_MPI_THREAD_MULTIPLE
};

extern int MPI_RANK;                                                            ///< Local MPI RANK in MPI parallel mode
extern int MPI_SIZE;                                                            ///< Total number of MPI RANKs in MPI parallel mode
extern int MPI_CART_RANK[3];                                                    ///< Cartesian coordinates of the current process if using MPI 3D domain decomposition
extern int MPI_CART_SIZE[3];                                                    ///< Number of processes used in each direction if using MPI 3D domain decomposition
extern bool MPI_3D_DECOMPOSITION;                                               ///< "true" if MPI should decompose in 3 dimensions

void* create_request();

void free_request(void* req);

int OP_MPI_Init_thread(int* argc, char ***argv, OP_MPI_Threads required, int* provided);

int OP_MPI_Comm_rank(OP_MPI_Comm commnicator, int* MPI_RANK);

int OP_MPI_Comm_size(OP_MPI_Comm commnicator, int* MPI_SIZE);

int OP_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                     OP_MPI_Datatype op_mpi_datatype, OP_MPI_Op op_mpi_op,
                     OP_MPI_Comm communicator);

int OP_MPI_Allreduce(OP_MPI_Options inplace, void *buf, int count,
                     OP_MPI_Datatype op_mpi_datatype, OP_MPI_Op op_mpi_op,
                     OP_MPI_Comm communicator);

int OP_MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                  OP_MPI_Datatype op_mpi_datatype, OP_MPI_Op op_mpi_op, int root,
                  OP_MPI_Comm communicator);

int OP_MPI_Isend(const void *buf, int count,  OP_MPI_Datatype op_mpi_datatype, int dest, int tag,
                 OP_MPI_Comm communicator, void *request);

int OP_MPI_Irecv(void *buf, int count,  OP_MPI_Datatype op_mpi_datatype, int dest, int tag,
                 OP_MPI_Comm communicator, void *request);

int OP_MPI_Wait(void *request, OP_MPI_Status status);

int OP_MPI_Cart_Setup(const int dims[], int& rank, int* cart_rank);

int OP_MPI_Bcast(void *buffer, int count,
                 OP_MPI_Datatype op_mpi_datatype, int root,
                 OP_MPI_Comm communicator);

void OP_MPI_Barrier(OP_MPI_Comm communicator);

void OP_MPI_Finalize();

void op_fftw_mpi_init();

fftw_plan op_fftw_mpi_plan_dft_r2c_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,
                                      double *in, fftw_complex *out,
                                      OP_MPI_Comm communicator,
                                      unsigned int flags);

fftw_plan op_fftw_mpi_plan_dft_c2r_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,
                                      fftw_complex *in, double *out,
                                      OP_MPI_Comm communicator,
                                      unsigned int flags);

ptrdiff_t op_fftw_mpi_local_size_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,
                                    OP_MPI_Comm communicator,
                                    ptrdiff_t *local_n0, ptrdiff_t *local_0_start);

void op_fftw_mpi_cleanup();

void op_mpi_write_vtk(const std::string Filename, std::stringstream& buffer,
    std::stringstream& hbuffer,
    std::stringstream& tbuffer);
    
void op_mpi_write_data(const std::string Filename, std::vector<std::string> buffer, int Size_X, int Size_Y, int Size_Z);

[[noreturn]] void abortInvalidMPI();

#endif
