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

class ReactorODEs : public FuncEval
{
public:
    /**
     * Constructor
     * @param[in] sol Solution object specifying initial system state.
     */
    ReactorODEs(shared_ptr<Solution> sol)
    {
        /* ---------------------- INITIALIZE MEMBER VARS ---------------------- */

        // pointer to the system's ThermoPhase object. updated by the solver during
        // simulation to provide iteration-specific thermodynamic properties.
        m_gas = sol->thermo();

        // pointer to the kinetics manager. provides iteration-specific species
        // production rates based on the current state of the ThermoPhase.
        m_kinetics = sol->kinetics();

        // the system's constant pressure, taken from the provided initial state.
        m_pressure = m_gas->pressure();

        // number of chemical species in the system.
        m_nSpecies = m_gas->nSpecies();

        // resize the vector<double> storage containers for species partial molar enthalpies
        // and net production rates. internal values are updated and used by the solver
        // per iteration.
        m_hbar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);

        // number of equations in the ODE system. a conservation equation for each
        // species, plus a single energy conservation equation for the system.
        m_nEqs = m_nSpecies + 1;
    }

    /**
     * Evaluate the ODE right-hand-side function, ydot = f(t,y).
     *   - overridden from FuncEval, called by the integrator during simulation.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] p sensitivity parameter vector, length nparams()
     *   - note: sensitivity analysis isn't implemented in this example
     */
    void eval(double t, double *y, double *ydot, double *p) override
    {
        // the solution vector *y* is [T, Y_1, Y_2, ... Y_K], where T is the
        // system temperature, and Y_k is the mass fraction of species k.
        // similarly, the time derivative of the solution vector, *ydot*, is
        // [dT/dt, Y_1/dt, Y_2/dt, ... Y_K/dt].
        // the following variables are defined for clear and convenient access
        // to these vectors:
        double temperature = y[0];
        double *massFracs = &y[1];
        double *dTdt = &ydot[0];
        double *dYdt = &ydot[1];

        /* ------------------------- UPDATE GAS STATE ------------------------- */
        // the state of the ThermoPhase is updated to reflect the current solution
        // vector, which was calculated by the integrator.
        m_gas->setMassFractions_NoNorm(massFracs);
        m_gas->setState_TP(temperature, m_pressure);

        /* ----------------------- GET REQ'D PROPERTIES ----------------------- */
        double rho = m_gas->density();
        double cp = m_gas->cp_mass();
        m_gas->getPartialMolarEnthalpies(&m_hbar[0]);
        m_kinetics->getNetProductionRates(&m_wdot[0]);

        /* -------------------------- ENERGY EQUATION ------------------------- */
        // the rate of change of the system temperature is found using the energy
        // equation for a closed-system constant pressure ideal gas:
        //     m*cp*dT/dt = - sum[h(k) * dm(k)/dt]
        // or equivalently:
        //     dT/dt = - sum[hbar(k) * dw(k)/dt] / (rho * cp)
        double hdot_vol = 0;
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            hdot_vol += m_hbar[k] * m_wdot[k];
        }
        *dTdt = -hdot_vol / (rho * cp);

        /* --------------------- SPECIES CONSERVATION EQS --------------------- */
        // the rate of change of each species' mass fraction is found using the closed-system
        // species conservation equation, applied once for each species:
        //     m*dY(k)/dt = dm(k)/dt
        // or equivalently:
        //     dY(k)/dt = dw(k)/dt * MW(k) / rho
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
        }
    }

    /**
     * Number of equations in the ODE system.
     *   - overridden from FuncEval, called by the integrator during initialization.
     */
    size_t neq() const override
    {
        return m_nEqs;
    }

    /**
     * Provide the current values of the state vector, *y*.
     *   - overridden from FuncEval, called by the integrator during initialization.
     * @param[out] y solution vector, length neq()
     */
    void getState(double *y) override
    {
        // the solution vector *y* is [T, Y_1, Y_2, ... Y_K], where T is the
        // system temperature, and Y_k is the mass fraction of species k.
        y[0] = m_gas->temperature();
        m_gas->getMassFractions(&y[1]);
    }

private:
    // private member variables, to be used internally.
    shared_ptr<ThermoPhase> m_gas;
    shared_ptr<Kinetics> m_kinetics;
    vector<double> m_hbar;
    vector<double> m_wdot;
    double m_pressure;
    size_t m_nSpecies;
    size_t m_nEqs;
};

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
    }
}

// The actual code is put into a function that
// can be called from the main program.
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
    fromxhat(x, ptcl, nx, rusr);

    //transfer arrays
    T[0] = ptcl[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y[ii - 1] = ptcl[ii];
    }

    //which Reactor ODE class?
    auto sol = newSolution("nDodecane_Reitz.yaml", "nDodecane_IG", "none");

    //set state
    auto gas = sol->thermo();
    gas->setState_TPY(T[0], p, Y);

    //set reactor and mechanisms
    auto odes = newReactor("Reactor", sol); 
    auto reactor = std::static_pointer_cast<Reactor>(odes); 
    ReactorNet net;
    net.addReactor(*reactor);

    //time and init
    double tnow = 0.0;
    net.setInitialTime(tnow);
    net.initialize();

    //integrate
    net.advance(dt);

    //stuff stuff
    size_t neq = reactor->neq();
    std::vector<double> y(neq);
    size_t n_species = gas->nSpecies();
    
    //get state and normalize
    reactor->getState(y.data());
    toxhat(y.data(), f, nx, rusr);
    myfnn(nx, x, fnn);

    //get reduced state
    for (int ii = 0; ii < nx; ii++)
    {
        f[ii] = f[ii] - x[ii] - fnn[ii];
    }

    //JACOBIAN START
    if (need[1] == 1)
    {
        Eigen::SparseMatrix<double> jac_spar = reactor->finiteDifferenceJacobian();
        Eigen::MatrixXd jac = Eigen::MatrixXd(jac_spar);

        for (size_t k = 0; k < n_species; k++)
        {
            g[k] = jac(k+1, 0);
        }
    }
    // {

    //     double xp[nx];
    //     double xm[nx];
    //     double fp[nf];
    //     double fm[nf];

    //     for (int ii = 0; ii < nx; ii++)
    //     {

    //         for (int jj = 0; jj < nx; jj++)
    //         {
    //             xp[jj] = x[jj];
    //             xm[jj] = x[jj];
    //         }

    //         xp[ii] += dx;
    //         xm[ii] -= dx;

    //         tnow = 0.0;
    //         fromxhat(xp, ptcl, nx, rusr);
    //         T[0] = ptcl[0];
    //         for (int ii = 1; ii < nx; ii++)
    //         {
    //             Y[ii - 1] = ptcl[ii];
    //         }
    //         gas->setState_TPY(T[0], p, Y);
    //         integrator->initialize(tnow, odes);
    //         integrator->integrate(dt);
    //         solution = integrator->solution();
    //         toxhat(solution, fp, nx, rusr);
    //         myfnn(nx, xp, fnn);
    //         for (int ii = 0; ii < nx; ii++)
    //         {
    //             fp[ii] = fp[ii] - xp[ii] - fnn[ii];
    //         }

    //         tnow = 0.0;
    //         fromxhat(xm, ptcl, nx, rusr);
    //         T[0] = ptcl[0];
    //         for (int ii = 1; ii < nx; ii++)
    //         {
    //             Y[ii - 1] = ptcl[ii];
    //         }
    //         gas->setState_TPY(T[0], p, Y);
    //         integrator->initialize(tnow, odes);
    //         integrator->integrate(dt);
    //         solution = integrator->solution();
    //         toxhat(solution, fm, nx, rusr);
    //         myfnn(nx, xm, fnn);
    //         for (int ii = 0; ii < nx; ii++)
    //         {
    //             fm[ii] = fm[ii] - xm[ii] - fnn[ii];
    //         }

    //         for (int jj = 0; jj < nx; jj++)
    //         {
    //             g[jj + ii * (nx)] = 1.0 * (fp[jj] - fm[jj]) / (2 * dx);
    //         }
    //     }
    // }
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

    d = alpha[0] * (T2 - T1);

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