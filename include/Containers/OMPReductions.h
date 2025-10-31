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

 *   File created :   2021
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef OMPREDUCTIONS_H
#define OMPREDUCTIONS_H

namespace openphase
{
#ifdef _OPENMP
    #pragma omp declare reduction (TensorD1Sum: Tensor<double, 1> : omp_out += omp_in) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (TensorD2Sum: Tensor<double, 2> : omp_out += omp_in) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (TensorD3Sum: Tensor<double, 3> : omp_out += omp_in) initializer (omp_priv = omp_orig)

    #pragma omp declare reduction (TensorV1Sum: Tensor<dVector3, 1> : omp_out += omp_in) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (TensorV2Sum: Tensor<dVector3, 2> : omp_out += omp_in) initializer (omp_priv = omp_orig)

    #pragma omp declare reduction (MatrixDMAX: Matrix<double> : omp_out = Matrix<double>::max(omp_in, omp_out)) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (MatrixDMIN: Matrix<double> : omp_out = Matrix<double>::min(omp_in, omp_out)) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (MatrixDSUM: Matrix<double> : omp_out += omp_in) initializer (omp_priv = omp_orig)

    #pragma omp declare reduction (vStressSUM: vStress : omp_out += omp_in) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (vStrainSUM: vStrain : omp_out += omp_in) initializer (omp_priv = omp_orig)

    #pragma omp declare reduction (dMatrix3x3SUM: dMatrix3x3 : omp_out += omp_in) initializer (omp_priv = omp_orig)
    #pragma omp declare reduction (dMatrix6x6SUM: dMatrix6x6 : omp_out += omp_in) initializer (omp_priv = omp_orig)

    #pragma omp declare reduction (dVectorNSUM: dVectorN : omp_out += omp_in) initializer (omp_priv = omp_orig)
#endif
}

#endif
