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
 *  Main contributors :   Reza Namdar
 *
 */
#include <numeric>
#include <stdexcept>
#include "ReactiveFlows/NeuralNetwork.h"

using namespace std;
using namespace openphase;

#ifdef ONNX

void NeuralNetwork::ReadInput(std::string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("NeuralNetwork");
    ConsoleOutput::WriteStandard("Source", InputFile);
    std::fstream inp(InputFile, std::ios::in);
    if (!inp) 
    {
        std::stringstream message;
        message << "File \"" << InputFile << "\" could not be opened";
        throw std::runtime_error(message.str());
    }
    std::stringstream inp_data; inp_data << inp.rdbuf(); inp.close();

    int moduleLocation = FileInterface::FindModuleLocation(inp_data, "NeuralNetwork");
    DO_PROP        = FileInterface::ReadParameterB(inp_data, moduleLocation, "DO_PROP",        false, false);
    DO_REACT       = FileInterface::ReadParameterB(inp_data, moduleLocation, "DO_REACT",       false, false);

    PropInputSize  = FileInterface::ReadParameterI(inp_data, moduleLocation, "PropInputSize",  false, 0);
    PropOutputSize = FileInterface::ReadParameterI(inp_data, moduleLocation, "PropOutputSize", false, 0);
    ReactInputSize = FileInterface::ReadParameterI(inp_data, moduleLocation, "ReactInputSize", false, 0);
    ReactOutputSize= FileInterface::ReadParameterI(inp_data, moduleLocation, "ReactOutputSize",false, 0);

    // NEW: optional paths (empty => disabled)
    PropModelPath  = FileInterface::ReadParameterS(inp_data, moduleLocation, "PropModelPath",  false, std::string{});
    ReactModelPath = FileInterface::ReadParameterS(inp_data, moduleLocation, "ReactModelPath", false, std::string{});
}

static size_t Product(const std::vector<int64_t>& s) 
{
    size_t p = 1;
    for (auto v : s) p *= static_cast<size_t>(v);
    return p;
}

void NeuralNetwork::GetFirstInOutNames(Ort::Session& s, std::string& in, std::string& out) 
{
    Ort::AllocatorWithDefaultOptions alloc;
    auto in0  = s.GetInputNameAllocated(0, alloc);
    auto out0 = s.GetOutputNameAllocated(0, alloc);
    in  = in0.get();
    out = out0.get();
}

void NeuralNetwork::LoadModels() 
{
    opts_ = Ort::SessionOptions{};
    opts_.SetIntraOpNumThreads(1);
    opts_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

    if (DO_PROP) 
    {
        if (PropModelPath.empty()) throw std::runtime_error("PropModelPath is empty");
        prop_sess_ = std::make_unique<Ort::Session>(env_, PropModelPath.c_str(), opts_);
        GetFirstInOutNames(*prop_sess_, prop_in_, prop_out_);
    }
    if (DO_REACT) 
    {
        if (ReactModelPath.empty()) throw std::runtime_error("ReactModelPath is empty");
        react_sess_ = std::make_unique<Ort::Session>(env_, ReactModelPath.c_str(), opts_);
        GetFirstInOutNames(*react_sess_, react_in_, react_out_);
    }
}

std::vector<float> NeuralNetwork::PredictProp(const std::vector<float>& x) 
{
    if (!prop_sess_) throw std::runtime_error("Prop model not loaded");
    if ((int)x.size() != PropInputSize) throw std::runtime_error("PredictProp: input size mismatch");

    Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    std::vector<int64_t> shape = {1, PropInputSize};

    Ort::Value in_tensor = Ort::Value::CreateTensor<float>(
        mem, const_cast<float*>(x.data()), x.size(), shape.data(), shape.size()
    );

    const char* in_names[]  = { prop_in_.c_str() };
    const char* out_names[] = { prop_out_.c_str() };

    auto outs = prop_sess_->Run(Ort::RunOptions{nullptr}, in_names, &in_tensor, 1, out_names, 1);

    auto& out0 = outs[0];
    auto out_shape = out0.GetTensorTypeAndShapeInfo().GetShape();
    size_t n = Product(out_shape);

    float* p = out0.GetTensorMutableData<float>();
    return std::vector<float>(p, p + n);
}

std::vector<float> NeuralNetwork::PredictReact(const std::vector<float>& x) 
{
    if (!react_sess_) throw std::runtime_error("React model not loaded");
    if ((int)x.size() != ReactInputSize) throw std::runtime_error("PredictReact: input size mismatch");

    Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    std::vector<int64_t> shape = {1, ReactInputSize};

    Ort::Value in_tensor = Ort::Value::CreateTensor<float>(
        mem, const_cast<float*>(x.data()), x.size(), shape.data(), shape.size()
    );

    const char* in_names[]  = { react_in_.c_str() };
    const char* out_names[] = { react_out_.c_str() };

    auto outs = react_sess_->Run(Ort::RunOptions{nullptr}, in_names, &in_tensor, 1, out_names, 1);

    auto& out0 = outs[0];
    auto out_shape = out0.GetTensorTypeAndShapeInfo().GetShape();
    size_t n = Product(out_shape);

    float* p = out0.GetTensorMutableData<float>();
    return std::vector<float>(p, p + n);
}

#endif
