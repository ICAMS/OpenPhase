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

 *   File created :   2016
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef SPHERICALHARMONICS_H
#define SPHERICALHARMONICS_H

#include "Includes.h"

namespace openphase
{

class SphericalHarmonics                                                        ///< Collection of real spherical harmonics up to six's order
{
    double operator()(int theta, int phi, double x, double y, double z)
    {
        int index = theta*10 + phi;
        double myReturn = 0.0;

        switch(index)
        {
            case 00:
            {
                myReturn =  Y00(x,y,z);
                break;
            }
            case 10:
            {
                myReturn =  Y10(x,y,z);
                break;
            }
            case 11:
            {
                myReturn =  Y11(x,y,z);
                break;
            }
            case 20:
            {
                myReturn =  Y20(x,y,z);
                break;
            }
            case 21:
            {
                myReturn =  Y21(x,y,z);
                break;
            }
            case 22:
            {
                myReturn =  Y22(x,y,z);
                break;
            }
            default:
            {
                ConsoleOutput::WriteExit("Non existing spherical harmonic implementation is called!" , "SphericalHarmonics", "operator()");
                OP_Exit(EXIT_FAILURE);
            }
        }
        return myReturn;
    }
    static double Y00(double x, double y, double z)
    {
        return 0.5*sqrt(1.0/Pi);
    }
    static double Y10(double x, double y, double z)
    {
        return sqrt(3.0/(4.0*Pi))*z;
    }
    static double Y11(double x, double y, double z)
    {
        return sqrt(3.0/(4.0*Pi))*x;
    }
    static double Y20(double x, double y, double z)
    {
        return 0.25*sqrt(5.0/Pi)*(3.0*z*z-1.0);
    }
    static double Y21(double x, double y, double z)
    {
        return 0.5*sqrt(15.0/Pi)*z*x;
    }
    static double Y22(double x, double y, double z)
    {
        return 0.25*sqrt(15.0/Pi)*(x*x - y*y);
    }

 protected:
 private:
};
    typedef double (SphericalHarmonics::*sh)(int theta, int phi);
}// namespace openphase
#endif
