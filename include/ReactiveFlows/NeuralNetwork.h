/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *  File created :   2025
 *  Main contributors :   Oleg Shchyglo; Reza Namdar
 *
 */
#ifndef NEURALNETWORK_H_INCLUDED
#define NEURALNETWORK_H_INCLUDED
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>

#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "BoundaryConditions.h"

#include "Tools/TimeInfo.h"
#include "RunTimeControl.h"
#include <vector>
#include <memory>
#include <stdexcept>

#ifdef ONNX

#include <onnxruntime_cxx_api.h>


using namespace std;

namespace openphase {
class OP_EXPORTS NeuralNetwork {
public:
    NeuralNetwork(const std::string& inputFileName, Ort::Env& env) : env_(env)
    {
        ReadInput(inputFileName);   // you already parse PropModelPath/ReactModelPath here
        LoadModels();               // new: create sessions
    }

    void ReadInput(std::string InputFile);

    bool DO_PROP{false};
    bool DO_REACT{false};

    int PropInputSize{0};
    int PropOutputSize{0};
    int ReactInputSize{0};
    int ReactOutputSize{0};

    // NEW: model paths
    std::string PropModelPath;   // e.g. "models/props.onnx"
    std::string ReactModelPath;  // e.g. "models/react.onnx"

    std::vector<float> PredictProp(const std::vector<float>& x);
    std::vector<float> PredictReact(const std::vector<float>& x);

    private:
    Ort::Env& env_;
    Ort::SessionOptions opts_{nullptr};

    std::unique_ptr<Ort::Session> prop_sess_;
    std::unique_ptr<Ort::Session> react_sess_;

    std::string prop_in_, prop_out_;
    std::string react_in_, react_out_;

    void LoadModels();
    static void GetFirstInOutNames(Ort::Session& s, std::string& in, std::string& out);
};
}

#endif

#endif


// ... existing includes/namespace ...



