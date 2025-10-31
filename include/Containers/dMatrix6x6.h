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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Philipp Engels
 *
 */

#ifndef DMATRIX6X6_H
#define DMATRIX6X6_H

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dMatrix3x3.h"
#include "dVector3.h"

namespace openphase
{

class dMatrix6x6
{
 public:
    dMatrix6x6()
    {
        memset(storage, 0, 36*sizeof(double));
    };
    dMatrix6x6(const dMatrix6x6& rhs)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] = rhs.storage[i][j];
        }
    };
    double& operator()(const size_t i, const size_t j)
    {
        assert(i < 6 && "Access beyond storage range");
        assert(j < 6 && "Access beyond storage range");

        return storage[i][j];
    };
    const double& operator()(const size_t i, const size_t j) const
    {
        assert(i < 6 && "Access beyond storage range");
        assert(j < 6 && "Access beyond storage range");

        return storage[i][j];
    };
    dMatrix6x6& set_to_zero(void)
    {
        memset(storage, 0, 36*sizeof(double));
        return *this;
    };
    dMatrix6x6& set_to_unity(void)
    {
        set_to_zero();
        for(int i = 0; i < 6; i++)
        {
            storage[i][i] = 1.0;
        }
        return *this;
    };

    dMatrix6x6& set_to_value(double value)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] = value;
        }
        return *this;
    };
    void pack(std::vector<double>& buffer)
    {
        for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
        {
            buffer.push_back(storage[i][j]);
        }
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
        {
            storage[i][j] = buffer[it]; ++it;
        }
    }
    double norm(void) const
    {
        double tmp = 0.0;
        for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
        {
            tmp += storage[i][j]*storage[i][j];
        }
        return sqrt(tmp);
    };
    /*double determinant(void) const
    {
        double determinant = 0.0;

        for (int i = 0; i < 6; i++)
        {
            double line_product = 1.0;
            for (int j = 0; j < 6; j++)
            {
                line_product *= storage[(i+j)%6][j];
            }
            determinant += line_product;
        }

        for (int i = 0; i < 6; i++)
        {
            double line_product = 1.0;
            for (int j = 0; j < 6; j++)
            {
                line_product *= storage[(i-j+6)%6][j];
            }
            determinant -= line_product;
        }
        std::cout << determinant << "\n";
        return determinant;
    };

    bool is_singular(void)
    {
        if(fabs(det()) > DBL_EPSILON)
        {
            return false;
        }
        else
        {
            return true;
        }
    };*/

    dMatrix6x6 operator*(const double m) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j]*m;
        }
        return tmp;
    };
    dMatrix6x6& operator*=(const double m)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] *= m;
        }
        return *this;
    };
    dMatrix6x6 operator*(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        tmp.set_to_zero();
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        for(int k = 0; k < 6; k++)
        {
            tmp(i,j) += storage[i][k]*rhs(k,j);
        }
        return tmp;
    };

    dMatrix6x6 operator+(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j] + rhs(i,j);
        }
        return tmp;
    };
    dMatrix6x6 operator-(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j] - rhs(i,j);
        }
        return tmp;
    };
    dMatrix6x6& operator+=(const dMatrix6x6& rhs)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] += rhs(i,j);
        }
        return *this;
    };
    dMatrix6x6 operator/(const double m) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j]/m;
        }
        return tmp;
    };
    dMatrix6x6& operator/=(const double m)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] /= m;
        }
        return *this;
    };
    dMatrix6x6& operator-=(const dMatrix6x6& rhs)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] -= rhs(i,j);
        }
        return *this;
    };
    dMatrix6x6& operator=(const dMatrix6x6& rhs)
    {
        memmove(data(), rhs.data(), 36*sizeof(double));
        return *this;
    };
    bool operator==(dMatrix6x6 rhs) const
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            if(fabs(storage[i][j] - rhs(i,j)) > DBL_EPSILON)
            {
                return false;
            }
        }
        return true;
    }
    bool operator!=(dMatrix6x6 rhs) const
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            if(fabs(storage[i][j] - rhs(i,j)) > DBL_EPSILON)
            {
                return true;
            }
        }
        return false;
    }
    dMatrix6x6 Hadamard_product(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j]*rhs(i,j);
        }
        return tmp;
    };
    dMatrix6x6& invert(void)
    {
        double Out[6][6];

        int indxc[6];
        int indxr[6];
        int ipiv[6] = {0, 0, 0, 0, 0, 0};
        int icol = 0;
        int irow = 0;
        double pivinv;
        double dum;
        double Uni[6][6] = {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            Out[i][j] = storage[i][j];
        }

        for(int i = 0; i < 6; i++)
        {
            double big = 0.0;
            for(int j = 0; j < 6; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < 6; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out[j][k]) >= big)
                    {
                        big = fabs(Out[j][k]);
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix."
                              << this->print() << "Matrix is Singular 1!!!\n";
                    OP_Exit(EXIT_FAILURE);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    double temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < 6; l++)
                {
                    double temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix.\n"
                          << this->print() << "Matrix is Singular 2!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
            for(int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
            for(int ll = 0; ll < 6; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = 5; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < 6; k++)
            {
                double temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] = Out[i][j];
        }
        return *this;
    };

    dMatrix6x6 inverted(void) const
    {
        double Out[6][6];

        int indxc[6];
        int indxr[6];
        int ipiv[6] = {0, 0, 0, 0, 0, 0};
        int icol = 0;
        int irow = 0;
        double pivinv;
        double dum;
        double Uni[6][6] = {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            Out[i][j] = storage[i][j];
        }

        for(int i = 0; i < 6; i++)
        {
            double big = 0.0;
            for(int j = 0; j < 6; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < 6; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out[j][k]) >= big)
                    {
                        big = fabs(Out[j][k]);
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix.\n"
                              << this->print() << "Matrix is Singular 1!!!\n";
                    OP_Exit(EXIT_FAILURE);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    double temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < 6; l++)
                {
                    double temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix.\n"
                          << this->print() << "Matrix is Singular 2!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
            for(int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
            for(int ll = 0; ll < 6; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = 5; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < 6; k++)
            {
                double temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        dMatrix6x6 tmp;
        tmp.set_to_zero();
        memmove(tmp.data(), Out, 36*sizeof(double));
        return tmp;
    };

    dMatrix6x6& transpose(void)
    {
        for(int i = 0; i < 5; i++)
        for(int j = i+1; j < 6; j++)
        {
            double tmp = storage[i][j];
            storage[i][j] = storage[j][i];
            storage[j][i] = tmp;
        }
        return *this;
    };

    dMatrix6x6 transposed(void) const
    {
        dMatrix6x6 tmp;

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[j][i];
        }
        return tmp;
    };

    dMatrix6x6& rotate(const dMatrix3x3& RotationMatrix)
    {
        double Out[3][3][3][3];

        for(int m = 0; m < 3; m++)
        for(int n = 0; n < 3; n++)
        for(int p = 0; p < 3; p++)
        for(int q = 0; q < 3; q++)
        {
            Out[m][n][p][q] = 0.0;

            for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
            for(int l = 0; l < 3; l++)
            {
                  // active rotation
                Out[m][n][p][q] += RotationMatrix(m,i)*
                                   RotationMatrix(n,j)*
                                   tensor(i,j,k,l)*
                                   RotationMatrix(p,k)*
                                   RotationMatrix(q,l);
            }
        }
        int VoigtIndex[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};
        for(int m = 0; m < 6; m++)
        for(int n = 0; n < 6; n++)
        {
            int i = VoigtIndex[m][0];
            int j = VoigtIndex[m][1];
            int k = VoigtIndex[n][0];
            int l = VoigtIndex[n][1];

            storage[m][n] = Out[i][j][k][l];
        }
        return *this;
    };
    dMatrix6x6 rotated(const dMatrix3x3& RotationMatrix)
    {
        double Out[3][3][3][3];

        for(int m = 0; m < 3; m++)
        for(int n = 0; n < 3; n++)
        for(int p = 0; p < 3; p++)
        for(int q = 0; q < 3; q++)
        {
            Out[m][n][p][q] = 0.0;

            for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
            for(int l = 0; l < 3; l++)
            {
                // active rotation
                Out[m][n][p][q] += RotationMatrix(m,i)*
                                   RotationMatrix(n,j)*
                                   tensor(i,j,k,l)*
                                   RotationMatrix(p,k)*
                                   RotationMatrix(q,l);
            }
        }
        dMatrix6x6 OUT;
        int VoigtIndex[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};

        for(int m = 0; m < 6; m++)
        for(int n = 0; n < 6; n++)
        {
            int i = VoigtIndex[m][0];
            int j = VoigtIndex[m][1];
            int k = VoigtIndex[n][0];
            int l = VoigtIndex[n][1];

            OUT(m,n) = Out[i][j][k][l];
        }
        return OUT;
    };
    dMatrix3x3 project(const dVector3& n)
    {
        dMatrix3x3 Out;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
        {
            Out(j,k) += n[i]*tensor(i,j,k,l)*n[l];
        }
        return Out;
    };
    std::string print(void) const
    {
        std::stringstream out;
        for(int i = 0; i < 6; i++)
        {
            out << "||" << std::setprecision(6)
                        << std::setw(10) << storage[i][0] << " "
                        << std::setw(10) << storage[i][1] << " "
                        << std::setw(10) << storage[i][2] << " "
                        << std::setw(10) << storage[i][3] << " "
                        << std::setw(10) << storage[i][4] << " "
                        << std::setw(10) << storage[i][5] << "||\n";
        }
        return out.str();
    };
    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 6; i++)
        {
           for(int j = 0; j < 6; j++)
           {
               out << storage[i][j] << sep;
           }
           out << "\n";
        }
        return out.str();
    };
    void read(std::fstream& inp)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            inp >> storage[i][j];
        }
    };
    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        for(int i = 0; i < 6; i++)
        {
           for(int j = 0; j < 6; j++)
           {
               out.push_back((double)storage[i][j]);
           }
        }
        return out;
    };
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        for(int i = 0; i < 6; i++)
        {
           for(int j = 0; j < 6; j++)
           {
               out.push_back((float)storage[i][j]);
           }
        }
        return out;
    };
    void read(std::stringstream& inp)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            inp >> storage[i][j];
        }
    };
    const double& tensor(const int i, const int j, const int k, const int l) const
    {
        return storage[(i==j)?(i):(6-(i+j))][(k==l)?(k):(6-(k+l))];
    };
    double* data(void)
    {
        return &storage[0][0];
    };
    const double* data(void) const
    {
        return &storage[0][0];
    };
    static dMatrix6x6 UnitTensor(void)
    {
        dMatrix6x6 myUnitTensor;
        return myUnitTensor.set_to_unity();
    };
    static dMatrix6x6 ZeroTensor(void)
    {
        dMatrix6x6 myZeroTensor;
        return myZeroTensor.set_to_zero();
    };
    constexpr size_t sizeX() const
    {
        return 6u;
    }
    constexpr size_t sizeY() const
    {
        return 6u;
    }
 protected:
 private:
    double storage[6][6];
};

}// namespace openphase
#endif
