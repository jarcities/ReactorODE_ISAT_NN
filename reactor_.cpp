

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
//CVODES stuff
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <cassert>
#include <iostream>
#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif

using namespace Cantera;
using namespace std;

namespace Gl
{

    shared_ptr<Solution> sol;
    shared_ptr<ThermoPhase> gas;

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

    }

}

using namespace Gl;

class ReactorODEs : public FuncEval
{
public:
    ReactorODEs(shared_ptr<Solution> sol)
    {
        m_gas = sol->thermo();
        m_kinetics = sol->kinetics();
        m_pressure = m_gas->pressure();
        m_nSpecies = m_gas->nSpecies();
        m_hbar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);
        m_nEqs = m_nSpecies + 1;
    }

    void eval(double t, double *y, double *ydot /*,double *p*/)
    // void eval(double t, double *y, double *ydot, double *p) override
    {
        double temperature = y[0];
        double *massFracs = &y[1];
        double *dTdt = &ydot[0];
        double *dYdt = &ydot[1];

        m_gas->setMassFractions_NoNorm(massFracs);
        m_gas->setState_TP(temperature, m_pressure);

        double rho = m_gas->density();
        double cp = m_gas->cp_mass();
        m_gas->getPartialMolarEnthalpies(&m_hbar[0]);
        m_kinetics->getNetProductionRates(&m_wdot[0]);

        //energy equation
        double hdot_vol = 0;
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            hdot_vol += m_hbar[k] * m_wdot[k];
        }
        *dTdt = -hdot_vol / (rho * cp);

        //species equation
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
        }
    }

    // //???
    // void eval(double t, const double *y, double *ydot)
    // {
    //     eval(t, const_cast<double*>(y), ydot, nullptr);
    // }

    size_t neq() const override
    {
        return m_nEqs;
    }

    void getState(double *y) override
    {
        y[0] = m_gas->temperature();
        m_gas->getMassFractions(&y[1]);
    }

    size_t nSpecies() const
    {
        return m_nSpecies;
    }

    const std::vector<std::string> &speciesNames() const
    {
        return m_gas->speciesNames();
    }

private:
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

    // for ( int ii = 0; ii < nx; ii++ ){ptcl[ii] = x[ii]*(100000.0*rusr[ii]);} 
}

void toxhat(double ptcl[], double x[], int &nx, double rusr[])
{

    x[0] = (ptcl[0] - rusr[0]) / rusr[nx];

    for (int ii = 1; ii < nx; ii++)
    {

        x[ii] = -log((ptcl[ii] + rusr[ii]) / rusr[ii]) / log(rusr[ii]);
    }

    // for ( int ii = 0; ii < nx; ii++ ){x[ii] = ptcl[ii]/(100000.0*rusr[ii]);}
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
        x1[ii] = x[ii];
    }

    for (int ll = 0; ll < nLayers; ll++)
    {

        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x2[kk] = 0.0;
            for (int jj = 0; jj < n1[ll]; jj++)
            {
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

static int RHS(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    ReactorODEs *odes = static_cast<ReactorODEs *>(user_data);
    double *y_data = NV_DATA_S(y);
    double *yd_data = NV_DATA_S(ydot);
    odes->eval((double)t, y_data, yd_data);
    return 0;
}

void myfgh(int need[], int &nx, double x[], int &nf, int &nh, int iusr[],
           double rusr[], double f[], double g[], double h[])
{

    double Y[nx - 1];
    double T[1];
    double ptcl[nx];
    // double *solution;
    double aTol = 1e-9;
    double rTol = 1e-6;
    double dt = rusr[2 * nx + 2];
    double dx = rusr[2 * nx + 3];
    double p = rusr[2 * nx + 4];
    // double fnn[nx];
    //set init conditions
    double tnow = 0.0;
    double t = tnow;
    const double tfinal = 1e-5;

    static int aaaa;

    if (aaaa != 7777)
    {
        Gl::initfgh();
        aaaa = 7777;
    }

    fromxhat(x, ptcl, nx, rusr);

    T[0] = ptcl[0];
    for (int ii = 1; ii < nx; ii++)
    {
        Y[ii - 1] = ptcl[ii];
    }

    gas->setState_TPY(T[0], p, Y);
    ReactorODEs odes = ReactorODEs(sol);
    size_t NEQ = odes.neq();
    // size_t Nsp = odes.nSpecies();
    // const auto &names = odes.speciesNames();

    //CVODE sens stuff
    SUNContext sunctx;
    int flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
    N_Vector y = N_VNew_Serial(NEQ, sunctx);
    assert(y);
    double *y_data = NV_DATA_S(y);
    odes.getState(y_data);
    void *m_cvode_mem = CVodeCreate(CV_BDF, sunctx); //line 287 cpp, line 98 header file
    assert(m_cvode_mem);
    flag = CVodeInit(m_cvode_mem, RHS, tnow, y);
    assert(flag >= 0);
    flag = CVodeSStolerances(m_cvode_mem, rTol, aTol);
    assert(flag >= 0);
    flag = CVodeSetUserData(m_cvode_mem, &odes);
    assert(flag >= 0);
    flag = CVodeSetMaxNumSteps(m_cvode_mem, 50000);
    assert(flag >= 0);
    flag = CVodeSetMaxStep(m_cvode_mem, dt); //1e-6 original
    assert(flag >= 0);
    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    assert(A);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    assert(LS);
    flag = CVodeSetLinearSolver(m_cvode_mem, LS, A);
    assert(flag >= 0);

    // shared_ptr<Integrator> integrator(newIntegrator("CVODE"));

    // integrator->initialize(tnow, odes);

    // integrator->setTolerances(aTol, rTol);

    // integrator->integrate(dt);

    // solution = integrator->solution();


    //JACOBIAN START
    N_Vector *yS = nullptr;
    // int Ns = 0;
    int Ns = (int)NEQ; 
    if (need[1] == 1)
    {
        //setup sensitivity analysis
        yS = N_VCloneVectorArray(Ns, y);
        assert(yS);
        for (int is = 0; is < Ns; ++is)
        {
            for (int j = 0; j < (int)NEQ; ++j)
            {
                /*init sens vectors as identity so each variables 
                init sens with respect to itself is 1 and for others 0*/
                NV_Ith_S(yS[is], j) = (is == j ? 1.0 : 0.0);
            }
        }
        //forward sens
        flag = CVodeSensInit(m_cvode_mem, Ns, CV_STAGGERED, /*fS*/ nullptr, yS); //line 235 cpp
        assert(flag >= 0);
        flag = CVodeSensEEtolerances(m_cvode_mem);
        assert(flag >= 0);
        vector<sunrealtype> p(NEQ);
        vector<sunrealtype> pbar(NEQ);
        for (int i = 0; i < NEQ; i++)
        {
            p[i] = y_data[i];                                            
            pbar[i] = (fabs(y_data[i]) > 1e-12) ? fabs(y_data[i]) : 1.0; 
        }
        //set sens parameters
        flag = CVodeSetSensParams(m_cvode_mem, p.data(), pbar.data(), nullptr);    
        assert(flag >= 0);
    }
    //JACOBIAN END

    //INTEGRATE
    flag = CVode(m_cvode_mem, tfinal, y, &t, CV_NORMAL);
    if (flag < 0)
    {
        std::cout << "integration failed -> " << flag << std::endl;
        if (yS) N_VDestroyVectorArray(yS, Ns);
        N_VDestroy(y);
        CVodeFree(&m_cvode_mem);
        SUNMatDestroy(A);
        SUNLinSolFree(LS);
        SUNContext_Free(&sunctx);
        return;
    }

    toxhat(NV_DATA_S(y), f, nx, rusr); 

    // myfnn(nx, x, fnn);

    for (int ii = 0; ii < nx; ii++)
    {
        // f[ii] = f[ii] - x[ii] - fnn[ii];
        f[ii] = f[ii] - x[ii];
    }

    //JACOBIAN START
    if (need[1] == 1)
    {
        //get final sens into jac matrix
        flag = CVodeGetSens(m_cvode_mem, &t, yS);
        assert(flag >= 0);

        for (int j = 0; j < (int)NEQ; ++j)
        {
            double *Sj = NV_DATA_S(yS[j]);
            for (int i = 0; i < (int)NEQ; ++i)
            {
                g[i + j * nx] = Sj[i];
            }
        }
    }
    //JACOBIAN END

    //deallocate
    if (yS) N_VDestroyVectorArray(yS, Ns);
    N_VDestroy(y);
    CVodeFree(&m_cvode_mem);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
    SUNContext_Free(&sunctx);
}

void mymix(int &nx, double x1[], double x2[], double alpha[], int iusr[], double rusr[])
{

    double Y1[nx - 1], Y2[nx - 1];
    double H1, H2;
    double T1[1], T2[1], d;
    double p = OneAtm;

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
    H2 -= alpha[0] * d;

    for (int ii = 0; ii < nx - 1; ii++)
    {
        d = Y2[ii] - Y1[ii];
        Y1[ii] += alpha[0] * d;
        Y2[ii] -= alpha[0] * d;
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