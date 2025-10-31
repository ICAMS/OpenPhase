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
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#ifndef GRAINSPROPERTIES_H
#define GRAINSPROPERTIES_H

#include "Includes.h"
#include "Settings.h"
#include "H5Interface.h"

namespace openphase
{

enum class GrainStages : int
{
    Seed,                                                                       ///< The state when the grain seed is just planted
    Nucleus,                                                                    ///< The state when the grain seed started to grow but not yet overcome the critical size
    Stable                                                                      ///< The state of a stable grain (above critical size)
};

enum class GrowthConstraintsViolations : int
{
    None,                                                                       ///< There is currently no violation
    Acute,                                                                      ///< An acute violation was detected in the last timestep, reduce mobility
    Recovering                                                                  ///< No violation was detected, but mobility is still recovering
};

class Grain                                                                     ///< Stores information about the grain (phase field). Used as the storage bit in the GrainInfo class
{
 public:
    bool Exist;                                                                 ///< True if the phase field exists in the simulation domain
    bool Mobile;                                                                ///< True if grain can move. Relevant for solid-liquid coupled motion
    size_t Parent;                                                              ///< Parent grain index
    size_t Phase;                                                               ///< Thermodynamic phase for a given phase field
    size_t Variant;                                                             ///< Symmetry variant of a given thermodynamic phase

    GrainStages  Stage;                                                         ///< Indicates growth stages of the grain (Seed, Nucleus, Stable)
    AggregateStates State;                                                      ///< Aggregate state of the grain (Solid, Liquid or Gas).

    double      Density;                                                        ///< Average density of the grain
    double      RefVolume;                                                      ///< Reference nucleus radius
    double      Volume;                                                         ///< Total volume of the grain in grid cells
    double      MAXVolume;                                                      ///< Maximum volume of the phase field during the simulation
    double      VolumeRatio;                                                    ///< 0.0 <= MAXVolume/RefVolume <= 1.0 ratio for supporting nuclei growth
    dVector3    Rcm;                                                            ///< Coordinates of the center of mass
    dVector3    Vcm;                                                            ///< Velocity of the center of mass
    dVector3    Acm;                                                            ///< Acceleration of the center of mass
    dVector3    aVel;                                                           ///< Angular velocity
    dVector3    aAcc;                                                           ///< Angular acceleration
    dVector3    Force;                                                          ///< Force acting on the grain
    dVector3    Torque;                                                         ///< Torque applied to the grain
    dMatrix3x3  InertiaM;                                                       ///< Tensor of the moment of inertia
    Quaternion  Orientation;                                                    ///< Orientation of the grain with respect to a reference frame

    GrowthConstraintsViolations GrowthConstraintsViolation;                     ///< Parameter used to detect violations of current grain growth constraints

    bool is_solid()  const {return State == AggregateStates::Solid;};
    bool is_liquid() const {return State == AggregateStates::Liquid;};
    bool is_gas()    const {return State == AggregateStates::Gas;};
    bool is_fluid()  const {return State != AggregateStates::Solid;};

    void Read(std::ifstream& inp)                                               ///< Reads grains info from a given file stream
    {
        inp.read(reinterpret_cast<char*>(&Exist        ), sizeof(bool));
        inp.read(reinterpret_cast<char*>(&Mobile       ), sizeof(bool));
        inp.read(reinterpret_cast<char*>(&Stage        ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Parent       ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Phase        ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Variant      ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&State        ), sizeof(int));
        inp.read(reinterpret_cast<char*>(&Density      ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&RefVolume    ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Volume       ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&MAXVolume    ), sizeof(double));
        VolumeRatio = MAXVolume/RefVolume;

        for (size_t i = 0; i < 3; i++)
        {
            inp.read(reinterpret_cast<char*>(&Rcm   [i]), sizeof(double));
            inp.read(reinterpret_cast<char*>(&Vcm   [i]), sizeof(double));
            inp.read(reinterpret_cast<char*>(&Acm   [i]), sizeof(double));
            inp.read(reinterpret_cast<char*>(&aVel  [i]), sizeof(double));
            inp.read(reinterpret_cast<char*>(&aAcc  [i]), sizeof(double));
            inp.read(reinterpret_cast<char*>(&Force [i]), sizeof(double));
            inp.read(reinterpret_cast<char*>(&Torque[i]), sizeof(double));

            for (size_t j = 0; j < 3; j++)
            {
                inp.read(reinterpret_cast<char*>(&InertiaM(i,j)), sizeof(double));
            }
        }

        for (size_t i = 0; i < 4; i++)
        {
            inp.read(reinterpret_cast<char*>(&Orientation[i]), sizeof(double));
        }
        Orientation.setRotationMatrix();
    }

    void ReadH5(std::vector<double> dbuffer, int& n)                            ///< Reads grains info from a given file stream
    {
        Exist         = dbuffer[n]; ++n;
        Mobile        = dbuffer[n]; ++n;
        Stage         = (GrainStages)(int)dbuffer[n]; ++n;
        Parent        = dbuffer[n]; ++n;
        Phase         = dbuffer[n]; ++n;
        Variant       = dbuffer[n]; ++n;
        State         = (AggregateStates)(int)dbuffer[n]; ++n;
        Density       = dbuffer[n]; ++n;
        RefVolume     = dbuffer[n]; ++n;
        Volume        = dbuffer[n]; ++n;
        MAXVolume     = dbuffer[n]; ++n;
        VolumeRatio   = MAXVolume/RefVolume;

        for (size_t i = 0; i < 3; i++)
        {
            Rcm   [i] = dbuffer[n]; ++n;
            Vcm   [i] = dbuffer[n]; ++n;
            Acm   [i] = dbuffer[n]; ++n;
            aVel  [i] = dbuffer[n]; ++n;
            aAcc  [i] = dbuffer[n]; ++n;
            Force [i] = dbuffer[n]; ++n;
            Torque[i] = dbuffer[n]; ++n;

            for (size_t j = 0; j < 3; j++)
            {
                InertiaM(i,j) = dbuffer[n]; ++n;
            }
        }

        for (size_t i = 0; i < 4; i++)
        {
            Orientation[i] = dbuffer[n]; ++n;
        }
        Orientation.setRotationMatrix();
    }

    void WriteH5(std::vector<double>& dbuffer)                                  ///< Reads grains info from a given file stream
    {
        dbuffer.push_back(Exist);
        dbuffer.push_back(Mobile);
        dbuffer.push_back((int)Stage);
        dbuffer.push_back(Parent);
        dbuffer.push_back(Phase);
        dbuffer.push_back(Variant);
        dbuffer.push_back((int)State);
        dbuffer.push_back(Density);
        dbuffer.push_back(RefVolume);
        dbuffer.push_back(Volume);
        dbuffer.push_back(MAXVolume);

        for (size_t i = 0; i < 3; i++)
        {
            dbuffer.push_back(Rcm   [i]);
            dbuffer.push_back(Vcm   [i]);
            dbuffer.push_back(Acm   [i]);
            dbuffer.push_back(aVel  [i]);
            dbuffer.push_back(aAcc  [i]);
            dbuffer.push_back(Force [i]);
            dbuffer.push_back(Torque[i]);

            for (size_t j = 0; j < 3; j++)
            {
                dbuffer.push_back(InertiaM(i,j));
            }
        }

        for (size_t i = 0; i < 4; i++)
        {
            dbuffer.push_back(Orientation[i]);
        }
    }

    void Write(std::ostream& outp) const                                        ///< Writes grains info to a given file stream
    {
        outp.write(reinterpret_cast<const char*>(&Exist        ), sizeof(bool));
        outp.write(reinterpret_cast<const char*>(&Mobile       ), sizeof(bool));
        outp.write(reinterpret_cast<const char*>(&Stage        ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Parent       ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Phase        ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Variant      ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&State        ), sizeof(int));
        outp.write(reinterpret_cast<const char*>(&Density      ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&RefVolume    ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Volume       ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&MAXVolume    ), sizeof(double));

        for (size_t i = 0; i < 3; i++)
        {
            outp.write(reinterpret_cast<const char*>(&Rcm   [i]), sizeof(double));
            outp.write(reinterpret_cast<const char*>(&Vcm   [i]), sizeof(double));
            outp.write(reinterpret_cast<const char*>(&Acm   [i]), sizeof(double));
            outp.write(reinterpret_cast<const char*>(&aVel  [i]), sizeof(double));
            outp.write(reinterpret_cast<const char*>(&aAcc  [i]), sizeof(double));
            outp.write(reinterpret_cast<const char*>(&Force [i]), sizeof(double));
            outp.write(reinterpret_cast<const char*>(&Torque[i]), sizeof(double));

            for (size_t j = 0; j < 3; j++)
            {
                outp.write(reinterpret_cast<const char*>(&InertiaM(i,j)), sizeof(double));
            }
        }

        for (size_t i = 0; i < 4; i++)
        {
            outp.write(reinterpret_cast<const char*>(&Orientation[i]), sizeof(double));
        }
    }

    void WriteTable(std::ostream& outp, const long int tStep, const char sep = ',',
                    const int precision = 16) const                             ///< Writes Statistics into csv Out
    {
        // Write header
        if (tStep == 0)
        {
            outp << "tStep"         << sep
                 << "Volume"        << sep
                 << "MAXVolume"     << sep
                 << "RefVolume"     << sep
                 << "Density"       << sep
                 << "Phase"         << sep
                 << "Parent"        << sep
                 << "X"         << sep << "Y"         << sep << "Z"         << sep
                 << "VX"        << sep << "VY"        << sep << "VZ"        << sep
                 << "AngularVX" << sep << "AngularVY" << sep << "AngularVZ" << sep
                 << "AX"        << sep << "AY"        << sep << "AZ"        << sep
                 << "AngularAZ" << sep << "AngularAY" << sep << "AngularAZ" << sep
                 << "ForceX"    << sep << "ForceY"    << sep << "ForceZ"    << sep
                 << "TorqueX"   << sep << "TorqueY"   << sep << "TorqueZ"   << sep
                 << "InertiaXX" << sep << "TorqueXY"  << sep << "TorqueXZ"  << sep
                 << "InertiaYX" << sep << "TorqueYY"  << sep << "TorqueYZ"  << sep
                 << "InertiaZX" << sep << "TorqueZY"  << sep << "TorqueZZ"  << sep
                 << "Q1"        << sep << "Q2"        << sep << "Q3"        << sep << "Q4\n";
        }

        outp << std::scientific  << std::setprecision(precision);
        outp << tStep         << sep;
        outp << Volume        << sep;
        outp << MAXVolume     << sep;
        outp << RefVolume     << sep;
        outp << Density       << sep;
        outp << Phase         << sep;
        outp << Parent        << sep;
        outp << Rcm   [0] << sep << Rcm   [1] << sep << Rcm   [2] << sep;
        outp << Vcm   [0] << sep << Vcm   [1] << sep << Vcm   [2] << sep;
        outp << Acm   [0] << sep << Acm   [1] << sep << Acm   [2] << sep;
        outp << aVel  [0] << sep << aVel  [1] << sep << aVel  [2] << sep;
        outp << aAcc  [0] << sep << aAcc  [1] << sep << aAcc  [2] << sep;
        outp << Force [0] << sep << Force [1] << sep << Force [2] << sep;
        outp << Torque[0] << sep << Torque[1] << sep << Torque[2] << sep;

        for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
        {
            outp << InertiaM(i,j) << sep;
        }

        for (size_t i = 0; i < 4; i++)
        {
            outp << Orientation[i] << sep;
        }
        outp << "\n";
    }

    void WriteTableVolumesHeader(std::ostream& outp, const char sep = ',') const ///< Writes Statistics into csv Out
    {
        outp << "Volume"    << sep
             << "MAXVolume" << sep
             << "RefVolume" << "\n";
    }

    void WriteTableVolumes(std::ostream& outp, const char sep = ',',
            const int precision = 16) const ///< Writes Statistics into csv Out
    {
        if (Volume != 0.0)
        {
            outp << std::scientific  << std::setprecision(precision);
            outp << Volume    << sep;
            outp << MAXVolume << sep;
            outp << RefVolume << "\n";
        }
    }

    Grain()                                                                     ///< Constructor - sets all values to their default values
    {
        Clear();
    }

    void Clear()                                                                ///< Clears the grains information
    {
        Exist         = false;
        Mobile        = true;
        Parent        = 0;
        Phase         = 0;
        Variant       = 0;
        Stage         = GrainStages::Stable;
        State         = AggregateStates::Solid;
        RefVolume     = 1.0;
        Volume        = 0.0;
        MAXVolume     = 0.0;
        Density       = 0.0;
        VolumeRatio   = 1.0;
        Orientation.set_to_zero();
        InertiaM.set_to_zero();
        Rcm.set_to_zero();
        Vcm.set_to_zero();
        Acm.set_to_zero();
        aVel.set_to_zero();
        aAcc.set_to_zero();
        Torque.set_to_zero();
        Force.set_to_zero();

        GrowthConstraintsViolation = GrowthConstraintsViolations::None;
    };

    void ClearForcesAndAccelerations()                                          ///< Clears the grains information
    {
        Acm.set_to_zero();
        aAcc.set_to_zero();
        Torque.set_to_zero();
        Force.set_to_zero();
    };

    std::vector<double>pack()                                                   ///< Creates the buffer for MPI communication
    {
        std::vector<double> loc_buffer(16,0.0);

        loc_buffer[0]  = (double)Exist;
        loc_buffer[1]  = (double)Mobile;
        loc_buffer[2]  = (double)Parent;
        loc_buffer[3]  = (double)Phase;
        loc_buffer[4]  = (double)Variant;
        loc_buffer[5]  = (double)Stage;
        loc_buffer[6]  = (double)State;
        loc_buffer[7]  =         RefVolume;
        loc_buffer[8]  =         Rcm[0];
        loc_buffer[9]  =         Rcm[1];
        loc_buffer[10] =         Rcm[2];
        loc_buffer[11] =         Orientation[0];
        loc_buffer[12] =         Orientation[1];
        loc_buffer[13] =         Orientation[2];
        loc_buffer[14] =         Orientation[3];

        loc_buffer[15] = (double)GrowthConstraintsViolation;

        return loc_buffer;
    }

    void unpack(std::vector<double>& loc_buffer)                                ///< Populates the Grain using the received MPI communication buffer
    {
        Exist          = (bool)  loc_buffer[0];
        Mobile         = (bool)  loc_buffer[1];
        Parent         = (size_t)loc_buffer[2];
        Phase          = (size_t)loc_buffer[3];
        Variant        = (size_t)loc_buffer[4];
        Stage          = (GrainStages)    loc_buffer[5];
        State          = (AggregateStates)loc_buffer[6];
        RefVolume      =         loc_buffer[7];
        Rcm[0]         =         loc_buffer[8];
        Rcm[1]         =         loc_buffer[9];
        Rcm[2]         =         loc_buffer[10];
        Orientation[0] =         loc_buffer[11];
        Orientation[1] =         loc_buffer[12];
        Orientation[2] =         loc_buffer[13];
        Orientation[3] =         loc_buffer[14];

        GrowthConstraintsViolation = (GrowthConstraintsViolations)loc_buffer[15];
    }

    bool IsNucleus(void) const
    {
        return Stage != GrainStages::Stable;
    };
    bool IsPresent(void) const
    {
        return Exist;
    };
};

class GrainsProperties
{
 public:
    GrainsProperties(void)
    {
    }
    GrainsProperties(const GrainsProperties& rhs)
    {
        GrainsStorage   = rhs.GrainsStorage;
    }
    void Allocate(const size_t size)
    {
        GrainsStorage.resize(size);
        for(size_t idx = 0; idx < GrainsStorage.size(); idx++)
        {
            GrainsStorage[idx].Clear();
        }
    }
    void Reallocate(const size_t size)
    {
        Resize(size);
    }
    void Resize(const size_t size)
    {
        size_t old_size = GrainsStorage.size();
        GrainsStorage.resize(size);
        for(size_t idx = old_size - 1; idx < GrainsStorage.size(); idx++)
        {
            GrainsStorage[idx].Clear();
        }
    }
    Grain& operator[](const size_t index)
    {
        assert (index < GrainsStorage.size());
        return GrainsStorage[index];
    }
    Grain const& operator[](const size_t index) const
    {
        assert (index < GrainsStorage.size());
        return GrainsStorage[index];
    }
    GrainsProperties& operator=(const GrainsProperties& rhs)
    {
        GrainsStorage   = rhs.GrainsStorage;
        return *this;
    }
    size_t size() const
    {
        return GrainsStorage.size();
    }
    double GrainRadius(const Grain& grain, int dim = 3) const
    {
        switch (dim)
        {
            case 1: return grain.Volume;
            case 2: return std::sqrt(grain.Volume/Pi);
            case 3: return std::pow(3.0/4.0*grain.Volume/Pi,1.0/3.0);
            default: return 0.0;
        }
    }
    size_t NumberOfGrains() const
    {
        size_t number = 0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0)
            {
                number++;
            }
        }
        return number;
    }
    size_t NumberOfGrains(size_t pIndex) const
    {
        size_t number = 0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                number++;
            }
        }
        return number;
    }
    double MeanGrainVolume(size_t pIndex) const
    {
        size_t number = 0;
        double volume = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume >= 1.0 && grain.Phase == pIndex)
            {
                volume+=grain.Volume;
                number++;
            }
        }
        if (number > 0) return volume/number;
        else return 0.0;
    }
    double MeanGrainRadius(size_t pIndex, int dim = 3) const
    {
        size_t number = 0;
        double mean = 0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                mean += GrainRadius(grain,dim);
                number++;
            }
        }
        if (number > 0) return mean/number;
        else return 0.0;
    }
    double VarianceOfGrainVolumes(size_t pIndex, double mu) const
    {
        size_t number = 0;
        double mu_2 = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                mu_2 += (grain.Volume-mu)*(grain.Volume-mu);
                number++;
            }
        }
        if (number > 0) return mu_2/number;
        else return 0.0;
    }
    double VarianceOfGrainVolumes(size_t pIndex) const
    {
        double mu = MeanGrainVolume(pIndex);
        return VarianceOfGrainVolumes(pIndex,mu);
    }
    double SkewnessOfGrainVolumes(size_t pIndex, double mu, double sigma) const
    {
        size_t number = 0;
        double mu_3 = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                mu_3 += (grain.Volume-mu)*(grain.Volume-mu)*(grain.Volume-mu);
                number++;
            }
        }
        if (number > 0 and sigma > 0) return mu_3/number/sigma/sigma/sigma;
        else return 0.0;
    }
    double SkewnessOfGrainVolumes(size_t pIndex) const
    {
        double mu = MeanGrainVolume(pIndex);
        double sigma = StandardDeviationOfGrainVolumes(pIndex);
        return SkewnessOfGrainVolumes(pIndex,mu,sigma);
    }
    template <typename T>
    double NthCentralMomentOfGrainVolumes(size_t pIndex, T N, double mu) const
    {
        size_t number = 0;
        double mu_N = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                mu_N += std::pow((grain.Volume-mu),N);
                number++;
            }
        }
        if (number > 0) return mu_N/number;
        else return 0.0;
    }
    template <typename T>
    double NthCentralMomentOfGrainVolumes(size_t pIndex, T N) const
    {
        double mu = MeanGrainVolume(pIndex);
        return NthCentralMomentOfGrainVolumes(pIndex,N,mu);
    }
    template <typename T>
    double NthStandardizedMomentOfGrainVolumes(size_t pIndex, T N, double mu, double sigma) const
    {
        return NthCentralMomentOfGrainVolumes(pIndex,N,mu)/std::pow(sigma,N);
    }
    template <typename T>
    double NthStandardizedMomentOfGrainVolumes(size_t pIndex, T N) const
    {
        double sigma = std::sqrt(VarianceOfGrainVolumes(pIndex));
        return NthCentralMomentOfGrainVolumes(pIndex,N)/std::pow(sigma,N);
    }
    double VarianceOfGrainRadii(size_t pIndex, double mu, int dim = 3) const
    {
        size_t number = 0;
        double mu_2 = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                double radius = GrainRadius(grain,dim);
                mu_2 += (radius-mu)*(radius-mu);
                number++;
            }
        }
        if (number > 0) return mu_2/number;
        else return 0.0;
    }
    double VarianceOfGrainRadii(size_t pIndex, int dim = 3) const
    {
        double mu = MeanGrainRadius(pIndex,dim);
        return VarianceOfGrainRadii(pIndex,mu,dim);
    }
    double SkewnessOfGrainRadii(size_t pIndex, double mu, double sigma, int dim = 3) const
    {
        size_t number = 0;
        double mu_3 = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                double radius = GrainRadius(grain,dim);
                mu_3 += (radius-mu)*(radius-mu)*(radius-mu);
                number++;
            }
        }
        if (number > 0) return mu_3/number/sigma/sigma/sigma;
        else return 0.0;
    }
    double SkewnessOfGrainRadii(size_t pIndex, int dim = 3) const
    {
        double mu = MeanGrainRadius(pIndex,dim);
        double sigma = StandardDeviationOfGrainRadii(pIndex,dim);
        return SkewnessOfGrainRadii(pIndex,mu,sigma,dim);
    }
    template <typename T>
    double NthCentralMomentOfGrainRadii(size_t pIndex, T N, double mu, int dim = 3) const
    {
        size_t number = 0;
        double mu_3 = 0.0;
        for (const Grain& grain : GrainsStorage)
        {
            //if (grain.Volume != 0.0) \\NOTE: counts invisible small grains
            if (grain.Volume > 1.0 && grain.Phase == pIndex)
            {
                double radius = GrainRadius(grain,dim);
                mu_3 += std::pow(radius-mu,N);
                number++;
            }
        }
        if (number > 0) return mu_3/number;
        else return 0.0;
    }
    template <typename T>
    double NthCentralMomentOfOfGrainRadii(size_t pIndex, T N, int dim = 3) const
    {
        double mu = MeanGrainRadius(pIndex,dim);
        return NthCentralMomentOfGrainRadii(pIndex,N,mu,dim);
    }
    template <typename T>
    double NthStandardizedMomentOfGrainRadii(size_t pIndex, T N, double mu, double sigma, int dim = 3) const
    {
        return NthCentralMomentOfGrainRadii(pIndex,N,mu,dim)/std::pow(sigma,N);
    }
    template <typename T>
    double NthStandardizedMomentOfGrainRadii(size_t pIndex, T N, int dim = 3) const
    {
        double sigma = StandardDeviationOfGrainRadii(pIndex,dim);
        return NthCentralMomentOfGrainRadii(pIndex,N,dim)/std::pow(sigma,N);
    }
    double StandardDeviationOfGrainVolumes(size_t pIndex) const
    {
        return std::sqrt(VarianceOfGrainVolumes(pIndex));
    }
    double StandardDeviationOfGrainRadii(size_t pIndex, bool dim = 3) const
    {
        return std::sqrt(VarianceOfGrainRadii(pIndex, dim));
    }
    double MeanGrainDiameter(size_t pIndex, bool dim = 3) const
    {
        return MeanGrainRadius(pIndex,dim)*2.0;
    }
    double StandardDeviationOfGrainDiameters(size_t pIndex, bool dim = 3) const
    {
        return StandardDeviationOfGrainRadii(pIndex,dim)*2.0;
    }
    double SkewnessOfGrainDiameters(size_t pIndex, bool dim = 3) const
    {
        return SkewnessOfGrainRadii(pIndex,dim)*2.0;
    }
    template <typename T>
    double NthCentralMomentOfGrainDiameters(size_t pIndex, T N, bool dim = 3) const
    {
        return NthCentralMomentOfOfGrainRadii(pIndex,N,dim)*2.0;
    }
    template <typename T>
    double NthStandardizedMomentOfGrainDiameters(size_t pIndex, T N, bool dim = 3) const
    {
        return NthStandardizedMomentOfGrainRadii(pIndex,N,dim)*2.0;
    }
    size_t add_grain(const size_t PhaseIndex)
    {
        for(size_t i = 0; i != GrainsStorage.size(); i++)
        if(not GrainsStorage[i].Exist)
        {
            GrainsStorage[i].Clear();
            GrainsStorage[i].Exist = true;
            GrainsStorage[i].Phase = PhaseIndex;
            GrainsStorage[i].Parent = 0;
            GrainsStorage[i].Orientation.set_to_zero();
            return i;
        }

        Grain locGrain;
        locGrain.Clear();
        locGrain.Exist  = true;
        locGrain.Phase  = PhaseIndex;
        locGrain.Parent = 0;
        locGrain.Orientation.set_to_zero();

        GrainsStorage.push_back(locGrain);
        return GrainsStorage.size() - 1;
    }
    bool PhasePresent(size_t pIndex) const
    {
        for(size_t idx = 0; idx < GrainsStorage.size(); idx++)
        {
            if(GrainsStorage[idx].Exist and
               GrainsStorage[idx].Phase == pIndex)
            {
                return true;
            }
        }
        return false;
    }
    void ChangePhaseIndex(size_t PhaseOld, size_t PhaseNew)
    {
        for(size_t idx = 0; idx < GrainsStorage.size(); idx++)
        if(GrainsStorage[idx].Phase == PhaseOld)
        {
            GrainsStorage[idx].Phase = PhaseNew;
        }
    }
    void ChangePhaseIndices(std::vector<size_t> PhaseOld, std::vector<size_t> PhaseNew)
    {
        if (PhaseOld.size() != PhaseNew.size())
        {
            std::stringstream message;
            message << "Incompatible phase indeces list sizes: "
                    << PhaseOld.size() << "!=" << PhaseNew.size() << "!" << "\n";
            ConsoleOutput::WriteExit(message.str(), thisclassname, "ChangePhaseIndices");
            OP_Exit(EXIT_FAILURE);
        }
        for(size_t idx = 0; idx < GrainsStorage.size(); idx++)
        for(size_t n = 0; n < PhaseOld.size();n++)
        if(GrainsStorage[idx].Phase == PhaseOld[n])
        {
            GrainsStorage[idx].Phase = PhaseNew[n];
        }
    }
    bool Write(const std::string FileName) const
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
        {
#endif
        std::ofstream outp(FileName.c_str(), std::ios::out | std::ios::binary);
        if (!outp)
        {
            ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be created!\n", thisclassname, "Write");
            return false;
        };

        const size_t size = GrainsStorage.size();
        outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

        for(size_t n = 0; n != size; n++)
        {
            GrainsStorage[n].Write(outp);
        }
        outp.close();
#ifdef MPI_PARALLEL
        }
#endif
        return true;
    }
    bool Write(const Settings& OPSettings, const int tStep) const
    {
        std::string FileName =
            FileInterface::MakeFileName(OPSettings.RawDataDir, thisclassname + "_", tStep,".dat");
        return Write(FileName);
    }
    bool Read(const std::string FileName)
    {
        std::ifstream inp(FileName.c_str(), std::ios::in | std::ios::binary);
        if (!inp)
        {
            ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
            return false;
        }

        size_t size = 0;
        inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));

        if(GrainsStorage.size() != size)
        {
            GrainsStorage.resize(size);
        }

        for(size_t n = 0; n != size; n++)
        {
            GrainsStorage[n].Read(inp);
        }

        inp.close();
        ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
        return true;
    }
    bool Read(const Settings& OPSettings, const int tStep)
    {
        std::string FileName = FileInterface::MakeFileName(OPSettings.InputRawDataDir, thisclassname + "_", tStep,".dat");
        return Read(FileName);
    }
    void WriteTable(const Settings& OPSettings, const size_t idx, const long int tStep, const char sep = ',',
            const int precision = 16) const
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
        {
#endif
        std::string FileName = OPSettings.TextDir + thisclassname + "_" + std::to_string(idx) + ".csv";
        std::ofstream outp(FileName, (tStep == 0) ? std::ios::app : std::ios::out);
        if (!outp)
        {
            ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be created!\n", thisclassname, "WriteTable");
            OP_Exit(EXIT_FAILURE);
        };

        GrainsStorage[idx].WriteTable(outp, tStep, sep, precision);
        outp.close();
#ifdef MPI_PARALLEL
        }
#endif
    }
    void WriteTable(const int tStep, const Settings& OPSettings, const char sep = ',',
            const int precision = 16) const
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
        {
#endif
        const std::string FileName = FileInterface::MakeFileName(OPSettings.TextDir, thisclassname + "_", tStep,".csv");
        std::ofstream outp(FileName);
        if (!outp)
        {
            ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be created!\n", thisclassname, "WriteTable");
            OP_Exit(EXIT_FAILURE);
        };

        for(const Grain& grain : GrainsStorage)
        {
            grain.WriteTable(outp, tStep, sep, precision);
        }
        outp.close();
#ifdef MPI_PARALLEL
        }
#endif
    }
    void WriteTableVolumes(const long int tStep, Settings& OPSettings,
            const char sep = ',', const int precision = 16) const
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
        {
#endif
        const std::string FileName = FileInterface::MakeFileName(OPSettings.TextDir, thisclassname + "_", tStep,".csv");
        std::ofstream outp(FileName);

        if (!outp)
        {
            std::stringstream message;
            message << "File \"" << FileName << "\" could not be opened\n";
            ConsoleOutput::WriteExit(message.str(), "GrainInfo", "WriteTableVolumes");
        };

        GrainsStorage[0].WriteTableVolumesHeader(outp, sep);
        for(const auto& grain : GrainsStorage)
        {
            grain.WriteTableVolumes(outp, sep, precision);
        }
        outp.close();
#ifdef MPI_PARALLEL
        }
#endif
    }
    //void Write(const size_t ID, const long int tStep) const
    //{
    //    std::ostringstream FileName;
    //    FileName.fill('0');
    //    FileName << "Checkpoint/PhaseField/GrainStats" << std::setw(8) << ID << "_" << tStep << ".dat";

    //    Write(FileName.str());
    //}
    //void Read(const size_t ID, const long int tStep)
    //{
    //    std::ostringstream FileName;
    //    FileName.fill('0');
    //    FileName << "Checkpoint/PhaseField/GrainStats" << std::setw(8) << ID << "_" << tStep << ".dat";

    //    Read(FileName.str());
    //}
    void WriteH5(const long int tStep, H5Interface& H5)
    {
        std::vector<double> data;
        data.push_back(GrainsStorage.size());
        for (size_t i = 0; i < GrainsStorage.size(); ++i)
        {
             GrainsStorage[i].WriteH5(data);
        }
        H5.WriteCheckPoint(tStep, thisclassname, data);
    }
    void ReadH5(const long int tStep, H5Interface& H5)
    {
        std::vector<double> data;
        H5.ReadCheckPoint(tStep, thisclassname, data);
        GrainsStorage.resize(data[0]);
        int n = 1;
        for (size_t i = 0; i < GrainsStorage.size(); ++i)
        {
             GrainsStorage[i].ReadH5(data,n);
        }
    }

    void Clear()
    {
        GrainsStorage.clear();
    }

    void ResetGrowthConstraintsViolations()
    {
        for(auto it = GrainsStorage.begin(); it != GrainsStorage.end(); ++it)
        {
            it->GrowthConstraintsViolation = GrowthConstraintsViolations::None;
        }
    }

    std::string thisclassname = "GrainsProperties";

 private:
    std::vector < Grain > GrainsStorage;
};

}// namespace openphase
#endif
