/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2020 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2020 OpenPhase Solutions GmbH,
 *                 Wasserstrasse 494, D-44795 Bochum, Germany.
 *
 *    All rights reserved.
 *
 *
 *    DEVELOPMENT VERSION, DO NOT PUBLISH OR DISTRIBUTE.
 *
 *
 *   OpenPhase (R) is a joint development of Interdisciplinary Centre for
 *   Advanced Materials Simulation (ICAMS), Ruhr University Bochum
 *   and OpenPhase Solutions GmbH.
 *
 *   File created :   2023
 *   Main contributors :   Marvin Tegeler; Oleg Shchyglo;
 *
 */

#include "mpi_wrapper.h"
#include "mpi.h"
#include "fftw3-mpi.h"
#include <string>
#include <cstdlib>
#include <signal.h>
#include <vector>
#include <sstream>
#include <iostream>

int MPI_RANK = 0;
int MPI_SIZE = 1;
int MPI_CART_RANK[3] = {0,0,0};
int MPI_CART_SIZE[3] = {1,1,1};
bool MPI_3D_DECOMPOSITION = false;

[[noreturn]] void abortInvalidMPI()
{
    std::cerr << "Invalid operation for MPI\n";
    std::abort();
}

MPI_Datatype getDatatype(OP_MPI_Datatype op_mpi_datatype)
{
    switch (op_mpi_datatype)
    {
    	case OP_MPI_Datatype::OP_MPI_CHAR: return MPI_CHAR;
        case OP_MPI_Datatype::OP_MPI_INT: return MPI_INT;
        case OP_MPI_Datatype::OP_MPI_LONG: return MPI_LONG;
        case OP_MPI_Datatype::OP_MPI_FLOAT: return MPI_FLOAT;
        case OP_MPI_Datatype::OP_MPI_DOUBLE: return MPI_DOUBLE;
        case OP_MPI_Datatype::OP_MPI_DOUBLE_INT: return MPI_DOUBLE_INT;
        case OP_MPI_Datatype::OP_MPI_UNSIGNED: return MPI_UNSIGNED;
        case OP_MPI_Datatype::OP_MPI_UNSIGNED_LONG: return MPI_UNSIGNED_LONG;
        case OP_MPI_Datatype::OP_MPI_UNSIGNED_LONG_LONG: return MPI_UNSIGNED_LONG_LONG;
        case OP_MPI_Datatype::OP_MPI_CXX_BOOL: return MPI_CXX_BOOL;
    }
    abortInvalidMPI();
}

MPI_Op getOp(OP_MPI_Op op_mpi_op)
{
    switch (op_mpi_op)
    {
        case OP_MPI_Op::OP_MPI_SUM: return MPI_SUM;
        case OP_MPI_Op::OP_MPI_PROD: return MPI_PROD;
        case OP_MPI_Op::OP_MPI_MIN: return MPI_MIN;
        case OP_MPI_Op::OP_MPI_MAX: return MPI_MAX;
        case OP_MPI_Op::OP_MPI_MAXLOC: return MPI_MAXLOC;
        case OP_MPI_Op::OP_MPI_LOR: return MPI_LOR;
        case OP_MPI_Op::OP_MPI_LAND: return MPI_LAND;
    }
    abortInvalidMPI();
}

void* create_request()
{
    MPI_Request* req = (MPI_Request*)malloc(sizeof(MPI_Request));
    return req;
}

void free_request(void* req)
{
    free(req);
}

int OP_MPI_Init_thread(int* argc, char ***argv, OP_MPI_Threads required, int* provided)
{
    int op_required = MPI_THREAD_SINGLE;

    switch(required)
    {
        case OP_MPI_THREAD_SINGLE:
            op_required = MPI_THREAD_SINGLE;
            break;
        case OP_MPI_THREAD_FUNNELED:
            op_required = MPI_THREAD_FUNNELED;
            break;
        case OP_MPI_THREAD_SERIALIZED:
            op_required = MPI_THREAD_SERIALIZED;
            break;
        case OP_MPI_THREAD_MULTIPLE:
            op_required = MPI_THREAD_MULTIPLE;
            break;
    }
    int result = MPI_Init_thread(argc, argv, op_required, provided);
    return result;
}

int OP_MPI_Comm_rank(OP_MPI_Comm communicator, int* MPI_RANK)
{
    int result = MPI_Comm_rank(MPI_COMM_WORLD, MPI_RANK);
    return result;
}

int OP_MPI_Comm_size(OP_MPI_Comm communicator, int* MPI_SIZE)
{
    int result = MPI_Comm_size(MPI_COMM_WORLD, MPI_SIZE);
    return result;
}

int OP_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                     OP_MPI_Datatype op_mpi_datatype, OP_MPI_Op op_mpi_op,
                     OP_MPI_Comm communicator)
{
    MPI_Datatype datatype = getDatatype(op_mpi_datatype);
    MPI_Op op = getOp(op_mpi_op);
    int result = MPI_Allreduce(sendbuf, recvbuf, count,
                               datatype, op, MPI_COMM_WORLD);
    return result;
}

int OP_MPI_Allreduce(OP_MPI_Options inplace, void *buf, int count,
                     OP_MPI_Datatype op_mpi_datatype, OP_MPI_Op op_mpi_op,
                     OP_MPI_Comm communicator)
{
    MPI_Datatype datatype = getDatatype(op_mpi_datatype);
    MPI_Op op = getOp(op_mpi_op);
    int result = MPI_Allreduce(MPI_IN_PLACE, buf, count,
                   datatype, op, MPI_COMM_WORLD);

    return result;
}

int OP_MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                  OP_MPI_Datatype op_mpi_datatype, OP_MPI_Op op_mpi_op, int root,
                  OP_MPI_Comm communicator)
{
    MPI_Datatype datatype = getDatatype(op_mpi_datatype);
    MPI_Op op = getOp(op_mpi_op);
    int result = MPI_Reduce(sendbuf, recvbuf, count,
                            datatype, op, root, MPI_COMM_WORLD);
    return result;
}

int OP_MPI_Isend(const void *buf, int count,
                 OP_MPI_Datatype op_mpi_datatype, int dest, int tag,
                 OP_MPI_Comm communicator, void *request)
{
    MPI_Datatype datatype = getDatatype(op_mpi_datatype);
    int result = MPI_Isend(buf, count, datatype, dest, tag,
                           MPI_COMM_WORLD, (MPI_Request*)request);
    return result;
}

int OP_MPI_Irecv(void *buf, int count,
                 OP_MPI_Datatype op_mpi_datatype, int dest, int tag,
                 OP_MPI_Comm communicator, void *request)
{
    MPI_Datatype datatype = getDatatype(op_mpi_datatype);
    int result = MPI_Irecv(buf, count, datatype, dest, tag,
                           MPI_COMM_WORLD, (MPI_Request*)request);
    return result;
}

int OP_MPI_Wait(void *request, OP_MPI_Status status)
{
    int result = MPI_Wait((MPI_Request*)request, MPI_STATUS_IGNORE);
    return result;
}

int OP_MPI_Cart_Setup(const int dims[], int& rank, int* cart_rank)
{
    MPI_Comm cart_comm;
    int reorder = 0;
    int periodic[3] = {0,0,0};

    int result1 = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periodic, reorder, &cart_comm);
    int result2 = MPI_Cart_coords(cart_comm, rank, 3, cart_rank);
    return std::max(result1,result2);
}

int OP_MPI_Bcast(void *buffer, int count,
                 OP_MPI_Datatype op_mpi_datatype, int root,
                 OP_MPI_Comm communicator)
{
    MPI_Datatype datatype = getDatatype(op_mpi_datatype);
    int result = MPI_Bcast(buffer, count, datatype, root, MPI_COMM_WORLD);
    return result;
}

void OP_MPI_Barrier(OP_MPI_Comm communicator)
{
    MPI_Barrier(MPI_COMM_WORLD);
}

void OP_MPI_Finalize()
{
    MPI_Finalize();
}

void op_fftw_mpi_init()
{
    fftw_mpi_init();
}

fftw_plan op_fftw_mpi_plan_dft_r2c_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,
                                      double *in, fftw_complex *out,
                                      OP_MPI_Comm communicator,
                                      unsigned int flags)
{
    return fftw_mpi_plan_dft_r2c_3d(n0, n1, n2, in, out, MPI_COMM_WORLD, flags);
}

fftw_plan op_fftw_mpi_plan_dft_c2r_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,
                                      fftw_complex *in, double *out,
                                      OP_MPI_Comm communicator,
                                      unsigned int flags)
{
    return fftw_mpi_plan_dft_c2r_3d(n0, n1, n2, in, out, MPI_COMM_WORLD, flags);
}

ptrdiff_t op_fftw_mpi_local_size_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2,
                                    OP_MPI_Comm communicator,
                                    ptrdiff_t *local_n0, ptrdiff_t *local_0_start)
{
    return fftw_mpi_local_size_3d(n0, n1, n2, MPI_COMM_WORLD, local_n0, local_0_start);
}

void op_fftw_mpi_cleanup()
{
    fftw_mpi_cleanup();
}

void op_mpi_write_vtk(const std::string Filename, std::stringstream& buffer,
    std::stringstream& hbuffer,
    std::stringstream& tbuffer)
{
    MPI_File fh;
    MPI_Status status;
    MPI_File_delete(Filename.c_str(), MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, Filename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset myoffset = 0;
    MPI_File_set_view(fh, myoffset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    size_t buffersize = buffer.str().size();
    size_t headersize = hbuffer.str().size();
    if (MPI_RANK == 0)
    {
        const std::string tmp = hbuffer.str();
        const char* headdata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, headdata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_Bcast(&headersize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    std::vector<size_t> sizes;
    sizes.resize(MPI_SIZE);
    MPI_Allgather(&buffersize, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    myoffset = headersize;
    for (int i = 0; i < MPI_RANK; ++i)
    {
        myoffset += sizes[i];
    }
    const std::string tmp = buffer.str();
    const char* bufferdata = tmp.c_str();
    MPI_File_write_at_all(fh, myoffset, bufferdata, tmp.size(), MPI_CHAR, &status);
    if (MPI_RANK == 0)
    {
        myoffset = headersize;
        for (int i = 0; i < MPI_SIZE; ++i)
        {
            myoffset += sizes[i];
        }
        const std::string tmp = tbuffer.str();
        const char* taildata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, taildata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_File_close(&fh);
}

void op_mpi_write_data(const std::string Filename, std::vector<std::string> buffer, int Size_X, int Size_Y, int Size_Z)
{
/*
    MPI_File fh;
    MPI_Status status;
    MPI_File_delete(Filename.c_str(), MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, Filename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset myoffset = 0;
    MPI_File_set_view(fh, myoffset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    size_t buffersize = buffer.size();
    std::vector<size_t> sizes(MPI_SIZE);
    MPI_Allgather(&buffersize, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    size_t maxsize= 0;
    for (int i = 0; i < MPI_SIZE; ++i)
    {
    	maxsize = std::max(maxsize, sizes[i]);
	}           	   
    std::vector<size_t> buffersizes(maxsize,0);    
    for (int i = 0; i < buffer.size(); ++i)
    {
    	buffersizes[i] =  buffer[i].size();
	}     
	std::vector<size_t> rbuffersizes(maxsize*MPI_SIZE,0);   
	std::vector<int> rsizes(MPI_SIZE, maxsize);           
	std::vector<int> displacements(MPI_SIZE);
	for (int i = 0; i < MPI_SIZE; ++i)
    {
    	displacements[i] = i*maxsize;
	}    
	MPI_Allgatherv(buffersizes.data(), maxsize , MPI_LONG_LONG,
                   rbuffersizes.data(), rsizes.data(), displacements.data(), MPI_LONG_LONG,
                  MPI_COMM_WORLD);
   	std::cout << Size_X << " " << Size_Y << " " << Size_Z << " " << MPI_CART_RANK[0] << " " << MPI_CART_RANK[1] << " " << MPI_CART_RANK[2] << " " << MPI_CART_SIZE[0] << " " << MPI_CART_SIZE[1] << " " << MPI_CART_SIZE[2] << std::endl; 	        
   	std::cout << "Displacements = " << std::endl;
   	for (int i = 0; i < displacements.size(); ++i)
   	{
  		std::cout << "d[" << i << "] = " << displacements[i] << std::endl;
  	}
  	std::cout << "rbuffersizes = " << std::endl;
   	for (int i = 0; i < rbuffersizes.size(); ++i)
   	{
  		std::cout << "r[" << i << "] = " << rbuffersizes[i] << std::endl;
  	}
    for (int i = 0; i < Size_Y; ++i)
    for (int j = 0; j < Size_Z; ++j)
    {
    	myoffset = 0;
    	for (int z = 0; z <= MPI_CART_RANK[2]; ++z)
    	for (int y = 0; y <= ((z<MPI_CART_RANK[2])?MPI_CART_SIZE[1]:MPI_CART_RANK[1]); ++y)
    	for (int x = 0; x <= ((z<MPI_CART_RANK[2] && y < MPI_CART_RANK[1])?MPI_CART_SIZE[0]:MPI_CART_RANK[0]); ++x)
		{
    		for (int r = 0; r <= ((z<MPI_CART_RANK[2] && y < MPI_CART_RANK[1] && x < MPI_CART_RANK[0])?Size_Z:j); ++r)
    		for (int q = 0; q <= ((z<MPI_CART_RANK[2] && y < MPI_CART_RANK[1] && x < MPI_CART_RANK[0] && r<j)?Size_Y:i); ++q)
    		{
		    	myoffset += rbuffersizes[displacements[x+MPI_CART_SIZE[0]*(y+MPI_CART_SIZE[1]*z)]+q+r*Size_Y];
		    	//std::cout << ((z<MPI_CART_RANK[2] && y < MPI_CART_RANK[1])?MPI_CART_SIZE[0]:MPI_CART_RANK[0]) << std::endl;
		    	//std::cout << q << " " << r << " " << x << " " << y << " " << z << " = " << displacements[x+MPI_CART_SIZE[0]*(y+MPI_CART_SIZE[1]*z)] << " " << rbuffersizes[displacements[x+MPI_CART_SIZE[0]*(y+MPI_CART_SIZE[1]*z)]+q+r*Size_Y] << " " << myoffset << std::endl;
		  	}
		}
		const std::string tmp = (i<buffer.size())?buffer[i]:"";
		const char* bufferdata = tmp.c_str();
		std::cout << tmp.size() << " " << myoffset << std::endl;
		MPI_File_write_at_all(fh, myoffset, bufferdata, tmp.size(), MPI_CHAR, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&fh);
    exit(0);*/
}
/*
void op_mpi_read_data(const std::string Filename, std::vector<int> window)
{
    MPI_File fh;
    MPI_Status status;
    MPI_File_delete(Filename.c_str(), MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, Filename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset myoffset = 0;
    MPI_File_set_view(fh, myoffset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    size_t buffersize = buffer.str().size();
    size_t headersize = hbuffer.str().size();
    if (MPI_RANK == 0)
    {
        const std::string tmp = hbuffer.str();
        const char* headdata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, headdata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_Bcast(&headersize, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    std::vector<size_t> sizes;
    sizes.resize(MPI_SIZE);
    MPI_Allgather(&buffersize, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    myoffset = headersize;
    for (int i = 0; i < MPI_RANK; ++i)
    {
        myoffset += sizes[i];
    }
    const std::string tmp = buffer.str();
    const char* bufferdata = tmp.c_str();
    MPI_File_write_at_all(fh, myoffset, bufferdata, tmp.size(), MPI_CHAR, &status);
    if (MPI_RANK == 0)
    {
        myoffset = headersize;
        for (int i = 0; i < MPI_SIZE; ++i)
        {
            myoffset += sizes[i];
        }
        const std::string tmp = tbuffer.str();
        const char* taildata = tmp.c_str();
        MPI_File_write_at(fh, myoffset, taildata, tmp.size(), MPI_CHAR, &status);
    }
    MPI_File_close(&fh);
}*/
