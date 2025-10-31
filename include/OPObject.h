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

 *   File created :   2014
 *   Main contributors :   Efim Borukhovich; Oleg Shchyglo
 *
 */

/*
 * OPObject is a parent virtual class for most of the stand-alone classes in
 * the OpenPhase. It defines the names of common methods and parameters which
 * unifies their syntax across the library. The common base also allows creating
 * lists of different objects and allows batch execution of common methods,
 * e.g. Remesh() and Advect() on a list of different objects.
 */

#ifndef OPOBJECT_H
#define OPOBJECT_H

#include "Includes.h"

namespace openphase
{

class Settings;
class BoundaryConditions;
class PhaseField;
class Velocities;
class AdvectionHR;

class OP_EXPORTS OPObject
{
 public:
    virtual ~OPObject(void);

    std::string thisclassname;                                                  ///< Object's implementation class name
    std::string thisobjectname;                                                 ///< Object's name

    bool initialized         = false;                                           ///< True if obect's Initialize has been executed
    bool input_read          = false;                                           ///< True if obect's ReadInput() has been executed
    bool remeshable          = false;                                           ///< True if the object has non-empty Remesh() method
    bool advectable          = false;                                           ///< True if the object has non-empty Advect() method
    bool readable            = false;                                           ///< True if the object has non-empty Read() method
    bool boundary_conditions = false;                                           ///< True if the object has non-empty SetBoundaryConditions()

    virtual void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "")///< Initializes internal variables and storages
    {
        (void) locSettings; //unused
    }

    virtual void ReadInput(const std::string InputFileName)                     ///< Reads input data from the user specified input file
    {
        (void) InputFileName; //unused
    }

    virtual void ReadInput(std::stringstream& inp)                              ///< Reads input data from the user specified stringstream
    {
        (void) inp; //unused
    }

    virtual void MoveFrame(const int dx, const int dy, const int dz,
                           const BoundaryConditions& BC)                        ///< Shifts the data in 3D storages by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction, correspondingly.
    {
            (void) dx; //unused
            (void) dy; //unused
            (void) dz; //unused
            (void) BC; //unused
    }

    virtual void Remesh(int newNx, int newNy, int newNz,
                        const BoundaryConditions& BC)                           ///< Remeshes the storages
    {
        (void) newNx; //unused
        (void) newNy; //unused
        (void) newNz; //unused
        (void) BC;    //unused
    }

    virtual void Advect(AdvectionHR& Adv, const Velocities& Vel,
                        PhaseField& Phi, const BoundaryConditions& BC,
                        const double dt, const double tStep)                    ///< Advects relevant fields
    {
        (void) Adv;   //unused
        (void) Vel;   //unused
        (void) Phi;   //unused
        (void) BC;    //unused
        (void) dt;    //unused
        (void) tStep; //unused
    }

    virtual bool Read(const Settings& locSettings,
                      const BoundaryConditions& BC, const int tStep)            ///< Reads binary input
    {
        (void) locSettings; //unused
        (void) BC;          //unused
        (void) tStep;       //unused
        return false;
    }

    virtual bool Write(const Settings& locSettings, const int tStep) const       ///< Writes binary output
    {
        (void) locSettings; //unused
        (void) tStep;       //unused
        return false;
    }

    virtual void SetBoundaryConditions(const BoundaryConditions& BC)
    {
        (void) BC; //unused
    }

    static OPObject* FindObject(std::vector<OPObject*> ObjectList,
                           std::string ObjectName, std::string thisclassname,
                           std::string thisfunctionname, bool necessary, bool verbose = true)
    /* Returns a reference to the first OPObject in the list whose class name
     * starts with ObjectName. I.e.: if ObjectName is "Elasticity", the first
     * object in the list which has either the name ElasticitySteinbach or
     * ElasticityKhachaturyan will be returned.
     */
    {
        for(unsigned int i = 0; i < ObjectList.size(); i++)
        {
            if(ObjectList[i]->thisclassname.size() >= ObjectName.size())
            {
                std::string locObjectName = ObjectList[i]->thisclassname;
                locObjectName.resize(ObjectName.size());

                if(locObjectName == ObjectName)
                {
                    return ObjectList[i];
                }
            }
        }
        if (verbose)
        {
            std::cout << "No " << ObjectName << " object found! " << thisclassname
                      << "." << thisfunctionname << "()" << std::endl;
        }

        if(necessary)
        {
            OP_Exit(13);
        }

        return nullptr;
    }

    static std::vector<OPObject*> FindObjects(std::vector<OPObject*> ObjectList,
                           std::string ObjectName, std::string thisclassname,
                           std::string thisfunctionname, bool necessary, bool verbose = true)
    /* Returns a list of OPObject references found in the passed ObjectList
     * whose class name start with ObjectName. I.e.: if ObjectName is
     * "Elasticity", all objects in the ObjectList which have a name starting
     * with "Elasticity", e.g. ElasticitySteinbach or ElasticityKhachaturyan
     * will be returned in the result list.
     */
    {
        std::vector<OPObject*> result;

        for(unsigned int i = 0; i < ObjectList.size(); i++)
        {
            if(ObjectList[i]->thisclassname.size() >= ObjectName.size())
            {
                std::string locObjectName = ObjectList[i]->thisclassname;
                locObjectName.resize(ObjectName.size());

                if(locObjectName == ObjectName)
                {
                    result.push_back(ObjectList[i]);
                }
            }
        }

        if(result.size() == 0)
        {
            if (verbose)
            {
                std::cout << "No " << ObjectName << " objects found! " << thisclassname
                          << thisfunctionname << std::endl;
            }
            if(necessary)
            {
                OP_Exit(13);
            }
        }
        return result;
    }

 protected:
 private:
};

}// namespace openphase
#endif
