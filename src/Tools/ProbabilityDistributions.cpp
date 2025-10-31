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

 *   File created :   2020
 *   Main contributors :   Oleg Shchyglo
 *
 */

#include "Tools/ProbabilityDistributions.h"

namespace openphase
{
using namespace std;

double ProbabilityDistributions::Normal(const double mu,
                                        const double sigma,
                                        const double value,
                                        const int derivative_order,
                                        const bool normalized)
{
    double norm = 1.0;
    if(normalized) norm = 1.0/(sigma*sqrt(2.0*Pi));
    double sigma_2 = 1.0/(sigma*sigma);
    double value1 = value - mu;
    double distribution = norm*exp(-0.5*pow(value1, 2)*sigma_2);
    double result = 0.0;

    switch(derivative_order)
    {
        case 0: // The distribution itself
        {
            result = distribution;
            break;
        }
        case 1: // First order derivative with respect to "value"
        {
            result = -distribution*value1*sigma_2;
            break;
        }
        case 2: // Second order derivative with respect to "value"
        {
            result = distribution*(pow(value1*sigma_2, 2) - sigma_2);
            break;
        }
        default:
        {
            //TODO: add warning and exit if non-existing distribution derivative is requested
        }
    }
    return result;
}

double ProbabilityDistributions::Cauchy(const double mu,
                                        const double gamma,
                                        const double value,
                                        const int derivative_order,
                                        const bool normalized)
{
    double norm = 1.0;
    if(normalized) norm = 1.0/(Pi*gamma);
    double gamma2 = gamma*gamma;
    double value1 = value - mu;
    double value2 = pow(value1, 2);
    double denominator_1 = 1.0/(value2 + gamma2);
    double distribution = norm*gamma2*denominator_1;
    double result = 0.0;

    switch(derivative_order)
    {
        case 0: // The distribution itself
        {
            result = distribution;
            break;
        }
        case 1: // First order derivative with respect to "value"
        {
            result = -2.0*norm*gamma2*value1*pow(denominator_1,2);
            break;
        }
        case 2: // Second order derivative with respect to "value"
        {
            result = 2.0*norm*gamma2*pow(denominator_1,2)*(4.0*value2*denominator_1 - 1.0);
            break;
        }
        default:
        {
            //TODO: add warning and exit if non-existing distribution derivative is requested
        }
    }
    return result;
}

}// namespace openphase
