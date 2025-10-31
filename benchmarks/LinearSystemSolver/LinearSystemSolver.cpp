#include "Includes.h"
#include "NumericalMethods/SystemOfLinearEquationsSolvers.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"

namespace op = openphase;

int main(int argc, char *argv[])
{
    std::string InputFileName = op::DefaultInputFileName;
    if (argc > 1) InputFileName = argv[1];

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    //feenableexcept(FE_OVERFLOW);   // the result was too large
    //feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::Settings OPSettings (InputFileName);
    op::BoundaryConditions  BC (OPSettings);
    op::TimeInfo            Timer;

    omp_set_num_threads(1);

    std::fstream inpF(InputFileName, std::ios::in);
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();
    const int  moduleLocation   = op::FileInterface::FindModuleLocation(inp, "LinearSystemSolver");
    const bool DoGauss          = op::FileInterface::ReadParameterB(inp, moduleLocation, std::string( "bGauss"));
    const bool DoJacobi         = op::FileInterface::ReadParameterB(inp, moduleLocation, std::string( "bJacobi"));
    const bool DoGaussSeidel    = op::FileInterface::ReadParameterB(inp, moduleLocation, std::string( "bGaussSeidel"));
    const bool DoGradientDecent = op::FileInterface::ReadParameterB(inp, moduleLocation, std::string( "bGradientDecent"));

    Timer.Initialize(OPSettings, "Execution Time Statistics");

    size_t Nx = OPSettings.Grid.Nx;
    size_t Ny = OPSettings.Grid.Ny;
    size_t Nz = OPSettings.Grid.Nz;

    size_t N = Nx*Ny*Nz;

    op::SparseMatrixCSR A(N,5);
    std::vector<double> b(N);
    std::vector<double> x(N);

    double dx = OPSettings.Grid.dx;

    auto idx = [Nx,Ny] (const int i, const int j, const int k) { return i+Nx*(j+Ny*k); };

    for (size_t i = 0; i < Nx; ++i)
    for (size_t j = 0; j < Ny; ++j)
    for (size_t k = 0; k < Nz; ++k)
    {
        x[idx(i,j,k)] = 0.00;
        double xx = (i-Nx/2.0)*dx;
        double yy = (j-Ny/2.0)*dx;
        double rr = sqrt(xx*xx + yy*yy);
        b[idx(i,j,k)] = 4.0*exp(-rr*rr)*(rr*rr -1);
        A(idx(i,j,k),idx(i,j,k)) = -4.0/dx/dx;

        const double Axp = 1.0/dx/dx;
        const double Axm = 1.0/dx/dx;
        const double Ayp = 1.0/dx/dx;
        const double Aym = 1.0/dx/dx;

        if (i+1 < Nx) A(idx(i,j,k),idx(i+1,j,k)) = Axp;
        //else switch (BC.BCNX)
        //{
        //    case op::NoFlux:   { A(idx(i,j,k),idx( i-1,j,k)) +=   Axp; break; }
        //    case op::Periodic: { A(idx(i,j,k),idx(   0,j,k)) +=   Axp; break; }
        //    case op::Fixed:    { A(idx(i,j,k),idx(   i,j,k)) +=   Axp; break; }
        //    default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        //}
        if (i >= 1) A(idx(i,j,k),idx(i-1,j,k)) = Axm;
        //else switch (BC.BC0X)
        //{
        //    case op::NoFlux:   { A(idx(i,j,k),idx( i+1,j,k)) +=   Axm; break; }
        //    case op::Periodic: { A(idx(i,j,k),idx(Nx-1,j,k)) +=   Axm; break; }
        //    case op::Fixed:    { A(idx(i,j,k),idx(   i,j,k)) +=   Axm; break; }
        //    default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        //}
        if (j+1 < Ny) A(idx(i,j,k),idx(i,j+1,k)) = Ayp;
        //else switch (BC.BCNY)
        //{
        //    case op::NoFlux:   { A(idx(i,j,k),idx(i,j-1,k)) +=   Ayp; break; }
        //    case op::Periodic: { A(idx(i,j,k),idx(i,  0,k)) +=   Ayp; break; }
        //    case op::Fixed:    { A(idx(i,j,k),idx(  i,j,k)) +=   Ayp; break; }
        //    default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        //}
        if (j >= 1) A(idx(i,j,k),idx(i,j-1,k)) = Aym;
        //else switch (BC.BC0Y)
        //{
        //    case op::NoFlux:   { A(idx(i,j,k),idx( i,j+1,k)) +=   Aym; break; }
        //    case op::Periodic: { A(idx(i,j,k),idx(i,Ny-1,k)) +=   Aym; break; }
        //    case op::Fixed:    { A(idx(i,j,k),idx(   i,j,k)) +=   Aym; break; }
        //    default: { throw std::invalid_argument("Boundary Condition not implemented!"); }
        //}
    }
    A.Sort();

    double MaxResidual = 1.0e-5;

    Timer.SetStart();

    op::ConsoleOutput::Write("Pivote Matrix");
    //op::SystemOfLinearEquationsSolvers::Pivote(A,b);
    op::SystemOfLinearEquationsSolvers::FastPivote(A,b);
    Timer.SetTimeStamp("Pivoting");

    double kappa = op::SystemOfLinearEquationsSolvers::SpectralNumberMAX(A,N);
    op::ConsoleOutput::Write("SpectralNumberMax", kappa);

    op::ConsoleOutput::Write("Precondition Jacobi");
    op::SystemOfLinearEquationsSolvers::PreconditionJacobi(A,b);
    Timer.SetTimeStamp("Precondition Jacobi");

    kappa = op::SystemOfLinearEquationsSolvers::SpectralNumberMAX(A,N);
    op::ConsoleOutput::Write("SpectralNumberMax", kappa);

    auto xGauss = x;
    if (DoGauss)
    {
        op::ConsoleOutput::Write("Calculate result with Gauss Method");
        auto Atmp = A; auto btmp = b;
        op::SystemOfLinearEquationsSolvers::Gauss(Atmp,xGauss,btmp);
        Timer.SetTimeStamp("Gauss Method");
    }

    auto xJacobi   = x;
    size_t iJacobi = 0;
    if (DoJacobi)
    {
        op::ConsoleOutput::Write("Calculate result with Jacobi's Method");
        iJacobi = op::SystemOfLinearEquationsSolvers::Jacobi(A,xJacobi,b,MaxResidual);
        Timer.SetTimeStamp("Jacobi Method");
    }

    auto xGaussSeidel   = x;
    size_t iGaussSeidel = 0;
    if (DoGaussSeidel)
    {
        op::ConsoleOutput::Write("Calculate result with Gauss-Seidel Method");
        iGaussSeidel = op::SystemOfLinearEquationsSolvers::GaussSeidel(A,xGaussSeidel,b,MaxResidual);
        Timer.SetTimeStamp("Gauss-Seidel Method");
    }

    auto xGradientDescent = x;
    size_t iGradientDescent = 0;
    if (DoGradientDecent)
    {
        op::ConsoleOutput::Write("Calculate result with Gradient Descent Method");
        op::SystemOfLinearEquationsSolvers::GradientDescent(A,xGradientDescent,b,MaxResidual);
        Timer.SetTimeStamp("Gradient Descent Method");
    }

    op::ConsoleOutput::Write("Calculate result with Conjugate Gradient Method");
    auto xConjugateGradient = x;
    size_t iConjugateGradient = op::SystemOfLinearEquationsSolvers::ConjugateGradient(A,xConjugateGradient,b,MaxResidual);
    Timer.SetTimeStamp("Conjugate Gradient Method");

    op::ConsoleOutput::Write("Calculate result with Biconjugate Gradient Method");
    auto xBiconjugateGradient = x;
    size_t iBiconjugateGradient = op::SystemOfLinearEquationsSolvers::BiconjugateGradient(A,xBiconjugateGradient,b,MaxResidual);
    Timer.SetTimeStamp("Biconjugate Gradient Method");

    //TODO write new BiCGStab solver the implemented fails!
    //auto xBiCGstab = x;
    //BiCGStab::BiCGstab_Jacobi(A,xBiCGstab,b,5,MaxResidual);
    //Timer.SetTimeStamp("xBiCGstab Method");

    Timer.PrintWallClockSummary();

    auto write_to_file = [&]<class T>(std::string name, T x)
    {
        const char sep = ',';
        double MaxDiff = 0.0;

        std::ofstream os(OPSettings.TextDir+'/'+name+".csv");
        os  << "i" << sep << name << sep << "Analytic" << "\n";
        for (size_t k = 0; k < Nz; ++k)
        for (size_t j = 0; j < Ny; ++j)
        for (size_t i = 0; i < Nx; ++i)
        {
            double xx = (i-Nx/2.0)*dx;
            double yy = (j-Ny/2.0)*dx;
            double rr = sqrt(xx*xx + yy*yy);
            double analytic = exp(-rr*rr);

            double numeric = x[idx(i,j,k)]  - x[idx(0,j,k)];

            os << i << sep << numeric << sep << analytic << "\n";

            if (std::abs(numeric  - analytic) > MaxDiff ) MaxDiff  = std::abs(numeric - analytic);
        }
        return MaxDiff;
    };

    const double DiffGauss               = write_to_file("Gauss"               , xGauss);
    const double DiffJacobi              = write_to_file("Jacobi"              , xJacobi);
    const double DiffGaussSeidel         = write_to_file("GaussSeidel"         , xGaussSeidel);
    const double DiffGradientDescent     = write_to_file("GradientDescent"      , xGradientDescent);
    const double DiffConjugateGradient   = write_to_file("ConjugateGradient"   , xConjugateGradient);
    const double DiffBiconjugateGradient = write_to_file("BiconjugateGradient" , xBiconjugateGradient);

    if (DoJacobi)         op::ConsoleOutput::Write("Iterations of Jacobi"               , iJacobi);
    if (DoGaussSeidel)    op::ConsoleOutput::Write("Iterations of Gauss-Seidel"         , iGaussSeidel);
    if (DoGradientDecent) op::ConsoleOutput::Write("Iterations of Gradient Descent"     , iGradientDescent);
    op::ConsoleOutput::Write("Iterations of Conjugate Gradient"   , iConjugateGradient);
    op::ConsoleOutput::Write("Iterations of Biconjugate Gradient" , iBiconjugateGradient);
    op::ConsoleOutput::Write("");
    if (DoGauss)          op::ConsoleOutput::Write("Max Error of Gauss"                , DiffGauss);
    if (DoJacobi)         op::ConsoleOutput::Write("Max Error of Jacobi"               , DiffJacobi);
    if (DoGaussSeidel)    op::ConsoleOutput::Write("Max Error of Gauss-Seidel"         , DiffGaussSeidel);
    if (DoGradientDecent) op::ConsoleOutput::Write("Max Error of Gradient Descent"     , DiffGradientDescent);
    op::ConsoleOutput::Write("Max Error of Conjugate Gradient"   , DiffConjugateGradient);
    op::ConsoleOutput::Write("Max Error of Biconjugate Gradient" , DiffBiconjugateGradient);

    std::ofstream os("Results.sim");
    os << std::scientific << std::setprecision(16);
    os << 3                          << "\n";
    if (DoGauss)          os << "DiffGauss "               << DiffGauss               << " " << DiffGauss               << "\n";
    if (DoJacobi)         os << "DiffJacobi "              << DiffJacobi              << " " << DiffJacobi              << "\n";
    if (DoGaussSeidel)    os << "DiffGaussSeidel "         << DiffGaussSeidel         << " " << DiffGaussSeidel         << "\n";
    if (DoGradientDecent) os << "DiffGradientDescent "     << DiffGradientDescent     << " " << DiffGradientDescent     << "\n";
    os << "DiffConjugateGradient "   << DiffConjugateGradient   << " " << DiffConjugateGradient   << "\n";
    os << "DiffBiconjugateGradient " << DiffBiconjugateGradient << " " << DiffBiconjugateGradient << "\n";
    os.close();

    return 0;
}
