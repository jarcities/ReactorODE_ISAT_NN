/*!
 * @file custom.cpp
 *
 * Custom Reactor
 *
 * Solve a closed-system constant pressure ignition problem where the governing
 * equations are custom-implemented, using Cantera's interface to CVODES to
 * integrate the equations.
 *
 * Keywords: combustion, reactor network, user-defined model, ignition delay
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "reactor.hpp"
#include <cmath>
#include "cantera/core.h"
#include "cantera/zerodim.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <chrono>
#include <string>
#include <iostream>
#include <cstdio>
#include <numeric>


using namespace Cantera;

namespace Gl
{

    shared_ptr<Solution> sol;    // = newSolution("h2o2.yaml", "ohmech", "none");
    shared_ptr<ThermoPhase> gas; // = sol->thermo();

    int nLayers = 5;
    int nNeurons = 10;
    int nx = 11;

    int ia[100];
    int ib[100];
    int n1[100];
    int n2[100];

    double A[1000000];
    double b[10000];

    void initfgh()
    {

        sol = newSolution("h2o2.yaml", "ohmech", "none");
        gas = sol->thermo();
    }

    void initfnn()
    {

        int i1 = 0;
        int i2 = 0;

        char file1[50];
        char file2[50];

        for (int ll = 1; ll <= nLayers; ll++)
        {

            ia[ll - 1] = i1;
            ib[ll - 1] = i2;

            sprintf(file1, "./A%d.csv", ll);
            sprintf(file2, "./B%d.csv", ll);

            n1[ll - 1] = nNeurons;
            n2[ll - 1] = nNeurons;
            if (ll == 1)
            {
                n1[ll - 1] = nx;
            }
            if (ll == nLayers)
            {
                n2[ll - 1] = nx;
            }

            FILE *pFile;
            float a;
            pFile = fopen(file1, "r+");
            for (int ii = 0; ii < n1[ll - 1] * n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                A[i1] = a;
                i1 = i1 + 1;
            }
            fclose(pFile);

            pFile = fopen(file2, "r+");
            for (int ii = 0; ii < n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                b[i2] = a;
                i2 = i2 + 1;
            }
            fclose(pFile);
        }

        // std::cerr<<A[ia[2]+3]<<std::endl;
        // std::cerr<<b[ib[1]+4]<<std::endl;
    }

}

using namespace Gl;

void fromxhat(double x[], double ptcl[], int &nx, double rusr[])
{

    ptcl[0] = (x[0] * rusr[nx]) + rusr[0];

    for (int ii = 1; ii < nx; ii++)
    {

        ptcl[ii] = rusr[ii] * exp(-log(rusr[ii]) * x[ii]) - rusr[ii];
    }

    // for ( int ii = 0; ii < nx; ii++ ){ptcl[ii] = x[ii]*(100000.0*rusr[ii]);} // REMOVE THIS
}

void toxhat(double ptcl[], double x[], int &nx, double rusr[])
{

    x[0] = (ptcl[0] - rusr[0]) / rusr[nx];

    for (int ii = 1; ii < nx; ii++)
    {

        x[ii] = -log((ptcl[ii] + rusr[ii]) / rusr[ii]) / log(rusr[ii]);
    }

    // for ( int ii = 0; ii < nx; ii++ ){x[ii] = ptcl[ii]/(100000.0*rusr[ii]);} //REMOVE THIS
}

double fAct(double x)
{

    return x * tanh(log(1.0 + exp(x)));
}

void myfnn(int &nx, double x[], double fnn[])
{

    static int bbbb;

    double x1[100];
    double x2[100];

    if (bbbb != 7777)
    {
        Gl::initfnn();
        bbbb = 7777;
    }

    for (int ii = 0; ii < n1[0]; ii++)
    {
        x1[ii] = x[ii]; // initialize the input
    }

    for (int ll = 0; ll < nLayers; ll++)
    {

        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x2[kk] = 0.0;
            for (int jj = 0; jj < n1[ll]; jj++)
            {
                // x2[kk] += A[ ia[ll] + jj + (kk-1)*n1[ll] ]*x1[jj];
                x2[kk] += A[ia[ll] + kk + (jj - 1) * n2[ll]] * x1[jj];
            }
            x2[kk] += b[ib[ll] + kk];

            if (ll < nLayers - 1)
            {
                x2[kk] = fAct(x2[kk]);
            }
        }

        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x1[kk] = x2[kk];
        }
    }

    for (int kk = 0; kk < nx; kk++)
    {
        fnn[kk] = x2[kk];
        // fnn[kk] = 0.0;
    }
}

void myfgh(int need[], int &nx, double x[], int &nf, int &nh, int iusr[],
           double rusr[], double f[], double g[], double h[])
{
    double Y[nx - 1]; //gas
    double T[1]; //temp
    double ptcl[nx]; 
    double *solution;
    double aTol = 1e-8; //rusr[2*nx];
    double rTol = 1e-8; //rusr[2*nx+1];
    double dt = rusr[2 * nx + 2];
    double dx = rusr[2 * nx + 3];
    double p = rusr[2 * nx + 4]; //pressure
    double fnn[nx];

    static int aaaa;
    if (aaaa != 7777)
    {
        Gl::initfgh();
        aaaa = 7777;
    }

    //un-normalize
    try
    {
        fromxhat(x, ptcl, nx, rusr);
    }
    catch (const std::exception& e)
    {
        std::cerr << "[myfgh] ERROR in fromxhat(): " << e.what() << "\n";
        std::cerr << "  x = [";
        for (int i = 0; i < nx; i++) {
            std::cerr << x[i] << (i+1<nx ? ", " : "");
        }
        std::cerr << "]\n";
        std::cerr << "  rusr = [";
        for (int i = 0; i < 2*nx+5; i++) {
            std::cerr << rusr[i] << (i+1<2*nx+5 ? ", " : "");
        }
        std::cerr << "]\n";
        throw;  
    }

    //transfer arrays
    T[0] = ptcl[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y[ii - 1] = ptcl[ii];
    }

    //set state
    try
    {
        gas->setState_TPY(T[0], p, Y);
    }
    catch (const Cantera::CanteraError& err) 
    {
        double sumY = std::accumulate(Y, Y + (nx-1), 0.0);
        std::cerr << "[myfgh] ERROR in setState_TPY(): " << err.what() << "\n";
        std::cerr << "  T = " << T[0] << ", p = " << p << "\n";
        std::cerr << "  Y = [";
        for (int i = 0; i < nx-1; i++) {
            std::cerr << Y[i] << (i+1<nx-1 ? ", " : "");
        }
        std::cerr << "], sum = " << sumY << "\n";
        throw;
    }

    //set reactor ODEs
    std::shared_ptr<Reactor> reactor;
    ReactorNet net;
    try
    {
        auto odes = newReactor("ConstPressureReactor", sol); 
        reactor = std::static_pointer_cast<Reactor>(odes); 
        net.addReactor(*reactor);
    }
    catch (const Cantera::CanteraError& err) 
    {
        std::cerr << "[myfgh] ERROR creating ConstPressureReactor: " 
                  << err.what() << "\n";
        std::cerr << "  sol ptr = " << sol.get() << "\n";
        throw;
    }

    //time and init
    try{
        double tnow = 0.0;
        net.setInitialTime(tnow);
        net.initialize();
    }
    catch (const Cantera::CanteraError& err)
    {
        std::cerr << "[myfgh] ERROR in net.initialize(): " << err.what() << "\n";
        std::cerr << "  initial time = 0.0\n";
        throw;
    }

    //integrate
    try 
    {
        net.advance(dt);
    } 
    catch (const Cantera::CanteraError& err) 
    {
        std::cerr << "!!! CVODE failed at dt="<<dt<<": "<<err.what()<<"\n";
        std::cerr << "    last x    = [";
        for (int i = 0; i < nx; i++) std::cerr << x[i] << (i+1<nx?", ":"");
        std::cerr << "]\n";
        std::cerr << "    last ptcl = [";
        for (int i = 0; i < nx; i++) std::cerr << ptcl[i] << (i+1<nx?", ":"");
        std::cerr << "]\n";
        throw;
    }

    //stuff stuff
    size_t neq = reactor->neq();
    std::vector<double> y(neq);
    size_t n_species = gas->nSpecies();
    
    //get state
    try
    {
        reactor->getState(y.data());
    }
    catch (const Cantera::CanteraError& err) 
    {
        std::cerr << "[myfgh] ERROR in getState(): " << err.what() << "\n";
        std::cerr << "  neq = " << reactor->neq() << "\n";
        throw;
    }
    
    //normalize and fnn
    try
    {
        toxhat(y.data(), f, nx, rusr);
        myfnn(nx, x, fnn);
    }
    catch (const std::exception& e)
    {
        std::cerr << "[myfgh] ERROR in toxhat()/myfnn(): " << e.what() << "\n";
        std::cerr << "  y  = [";
        for (int i = 0; i < nx; i++) std::cerr << y[i] << (i+1<nx?", ":"");
        std::cerr << "]\n";
        std::cerr << "  f  = [";
        for (int i = 0; i < nx; i++) std::cerr << f[i] << (i+1<nx?", ":"");
        std::cerr << "]\n";
        std::cerr << "  fnn= [";
        for (int i = 0; i < nx; i++) std::cerr << fnn[i] << (i+1<nx?", ":"");
        std::cerr << "]\n";
        throw;
    }

    //get reduced state
    for (int ii = 0; ii < nx; ii++)
    {
        f[ii] = f[ii] - x[ii] - fnn[ii];
    }

    //JACOBIAN START
    if (need[1] == 1)
    {
        Eigen::SparseMatrix<double> jac_sparse = reactor->finiteDifferenceJacobian();
        Eigen::MatrixXd jac = Eigen::MatrixXd(jac_sparse);

        Eigen::MatrixXd jac_nn(nx, nx);
        const double eps = dx;  
        std::vector<double> fnn_base(nx), fnn_p(nx), fnn_m(nx);
        myfnn(nx, x, fnn_base.data());

        for (int j = 0; j < nx; ++j) 
        {
            std::vector<double> x_p(x, x + nx), x_m(x, x + nx);
            x_p[j] += eps;
            x_m[j] -= eps;
            myfnn(nx, x_p.data(), fnn_p.data());
            myfnn(nx, x_m.data(), fnn_m.data());
            for (int i = 0; i < nx; ++i) {
                jac_nn(i, j) = (fnn_p[i] - fnn_m[i]) / (2 * eps);
                double jac_eig = jac(i, j);
                double identity = (i == j ? 1.0 : 0.0);
                g[i + j * nx] = jac_eig - identity - jac_nn(i, j);
            }
        }
    }
    //JACOBIAN END
}

void mymix(int &nx, double x1[], double x2[], double alpha[], int iusr[], double rusr[])
{

    double Y1[nx - 1], Y2[nx - 1];
    double H1, H2;
    double T1[1], T2[1], d;
    double p = OneAtm; // rusr[2*nx+4];

    T1[0] = x1[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y1[ii - 1] = x1[ii];
    }
    T2[0] = x2[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y2[ii - 1] = x2[ii];
    }

    gas->setState_TPY(T1[0], p, Y1);

    H1 = gas->enthalpy_mass();

    gas->setState_TPY(T2[0], p, Y2);

    H2 = gas->enthalpy_mass();

    d = H2 - H1;
    H1 += alpha[0] * d;
    H2 -= alpha[0] * d; // mix enthalpies

    for (int ii = 0; ii < nx - 1; ii++)
    {
        d = Y2[ii] - Y1[ii];
        Y1[ii] += alpha[0] * d;
        Y2[ii] -= alpha[0] * d; // mix mass fractions
    }

    // d = alpha[0] * (T2 - T1); //original
    d = alpha[0] * (T2[0] - T1[0]);

    gas->setState_TPY(T1[0] + d, p, Y1);
    gas->setState_HP(H1, p);

    T1[0] = gas->temperature();

    gas->setState_TPY(T2[0] - d, p, Y2);
    gas->setState_HP(H2, p);

    T2[0] = gas->temperature();

    x1[0] = T1[0];
    for (int ii = 1; ii < nx; ii++)
    {
        x1[ii] = Y1[ii - 1];
    }
    x2[0] = T2[0];
    for (int ii = 1; ii < nx; ii++)
    {
        x2[ii] = Y2[ii - 1];
    }
}