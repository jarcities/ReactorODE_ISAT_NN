#include "cantera/core.h"
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

#ifndef SUN_COMM_NULL
#define SUN_COMM_NULL NULL
#endif

namespace fs = std::filesystem;
using namespace Cantera;
using namespace std;

class ReactorODEs
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

    void eval(double t, const double *y, double *ydot /*, double *p */)
    {
        double temperature = y[0];
        const double *Y = y + 1;

        //update gas state
        m_gas->setMassFractions_NoNorm(Y);
        m_gas->setState_TP(temperature, m_pressure);

        //get properties
        double rho = m_gas->density();
        double cp = m_gas->cp_mass();
        m_gas->getPartialMolarEnthalpies(m_hbar.data());
        m_kinetics->getNetProductionRates(m_wdot.data());

        //energy equation
        double hdot_vol = 0.0;
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            hdot_vol += m_hbar[k] * m_wdot[k];
        }
        ydot[0] = -hdot_vol / (rho * cp);

        //species equation
        for (size_t k = 0; k < m_nSpecies; k++)
        {
            ydot[1 + k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
        }
    }

    size_t neq() const
    {
        return m_nEqs;
    }

    void getState(double *y)
    {
        y[0] = m_gas->temperature();
        m_gas->getMassFractions(y + 1);
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
    double m_pressure;
    size_t m_nSpecies, m_nEqs;
    std::vector<double> m_hbar, m_wdot;
};

static int RHS(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    ReactorODEs *odes = static_cast<ReactorODEs *>(user_data);
    double *y_data = NV_DATA_S(y);
    double *yd_data = NV_DATA_S(ydot);
    odes->eval((double)t, y_data, yd_data);
    return 0;
}

int main()
{
    /////////////////////////
    //// INIT CONDITIONS ////
    /////////////////////////
    std::string mechanism_ = "h2o2.yaml";
    std::string name_ = "ohmech";
    std::string comp_mix = "H2:2, O2:1, N2:4";
    // std::string mechanism_ = "nDodecane_Reitz.yaml";
    // std::string name_ = "nDodecane_IG";
    // std::string comp_mix = "c12h26:1, o2:1, n2:3.76";
    double temp_ = 1001.0;
    double pressure_ = OneAtm;
    const double t0 = 0.0;
    double t = t0;
    const double tfinal = 1e-4;
    const double reltol = 1e-6;
    const double abstol = 1e-9;
    int ToT = 1; //0 print, 1 save to csv, 2 both

    //init gas object and IC
    auto sol = newSolution(mechanism_, name_, "none");
    auto gas = sol->thermo();
    gas->setState_TPX(temp_, pressure_, comp_mix);

    //init reactor object
    ReactorODEs odes(sol);
    size_t NEQ = odes.neq();
    size_t Nsp = odes.nSpecies();
    const auto &names = odes.speciesNames();

    //////////////////////
    //// CVODES STUFF ////
    //////////////////////
    //REFERENCES
    //https://cantera.org/dev/cxx/d7/dd9/CVodesIntegrator_8h_source.html#l00098
    //https://cantera.org/dev/cxx/da/d4f/CVodesIntegrator_8cpp_source.html

    //create sundials object
    SUNContext sunctx;
    int flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);

    //allocate state (solution) vector
    N_Vector y = N_VNew_Serial(NEQ, sunctx);
    assert(y);
    double *y_data = NV_DATA_S(y);
    odes.getState(y_data);

    //create CVODES solver
    //"CV_BDF  - Use BDF methods"
    //"CV_NEWTON - use Newton's method"
    void *m_cvode_mem = CVodeCreate(CV_BDF, sunctx); //line 287 cpp, line 98 header file
    assert(m_cvode_mem);
    flag = CVodeInit(m_cvode_mem, RHS, t0, y);
    assert(flag >= 0);

    //set tolerances and step size
    flag = CVodeSStolerances(m_cvode_mem, reltol, abstol);
    assert(flag >= 0);
    flag = CVodeSetUserData(m_cvode_mem, &odes);
    assert(flag >= 0);
    flag = CVodeSetMaxNumSteps(m_cvode_mem, 50000);
    assert(flag >= 0);
    flag = CVodeSetMaxStep(m_cvode_mem, 1e-6);
    assert(flag >= 0);

    //set linear solver
    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    assert(A);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    assert(LS);
    flag = CVodeSetLinearSolver(m_cvode_mem, LS, A);
    assert(flag >= 0);

    //setup sensitivity analysis
    int Ns = (int)NEQ; 
    N_Vector *yS = N_VCloneVectorArray(Ns, y);
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
    flag = CVodeSetSensParams(m_cvode_mem,
                              p.data(),    
                              pbar.data(), 
                              nullptr);    
    assert(flag >= 0);

    //INTEGRATE
    flag = CVode(m_cvode_mem, tfinal, y, &t, CV_NORMAL);
    if (flag < 0)
    {
        cout << "integration failed -> " << flag << endl;
        return 1;
    }

    //get final sens into jac matrix
    flag = CVodeGetSens(m_cvode_mem, &t, yS);
    assert(flag >= 0);
    std::vector<std::vector<double>> J;
    J.resize(NEQ, std::vector<double>(NEQ));
    for (int j = 0; j < (int)NEQ; ++j)
    {
        double *Sj = NV_DATA_S(yS[j]);
        for (int i = 0; i < (int)NEQ; ++i)
        {
            J[i][j] = Sj[i];
        }
    }

    /////////////////
    //// RESULTS ////
    /////////////////
    if (ToT == 0 || ToT == 2)
    {
        //print init conditions and sens params
        cout << "init conditions:" << endl;
        cout << "==================" << endl;
        cout << "T = " << temp_ << " K" << endl;
        for (size_t k = 0; k < Nsp; k++)
        {
            double Yk = gas->massFraction(k); //get mass fraction
            cout << "Y_" << names[k] << ": " << Yk << endl;
        }
        cout << "no. of equations = " << NEQ << endl;
        cout << "no. of sens params = " << Ns << endl;

        //print final sol
        cout << "\nsolution:" << endl;
        cout << "===========" << endl;
        cout << "T = " << y_data[0] << " K" << endl;
        for (size_t k = 0; k < Nsp; k++)
        {
            cout << "Y_" << names[k] << ": " << y_data[k + 1] << endl;
        }
        
        //print jacobian
        cout << "\njacobian df[t]/dy[0]:" << endl;
        cout << "======================" << endl;
        cout << setw(12) << " ";
        cout << setw(13) << "dT_0";
        for (size_t k = 0; k < Nsp; k++)
        {
            cout << setw(15) << ("dY" + names[k] + "_0");
        }
        cout << endl;
        for (size_t i = 0; i < NEQ; i++)
        {
            if (i == 0)
            {
                cout << setw(10) << "dT_f:";
            }
            else
            {
                cout << setw(10) << ("dY" + names[i - 1] + "_f:");
            }
            for (size_t j = 0; j < NEQ; j++)
            {
                cout << setw(15) << scientific << setprecision(6) << J[i][j];
            }
            cout << endl;
        }
        //diag (self-sensitivities)
        cout << "\ndiag elements:" << endl;
        cout << "================" << endl;
        for (size_t i = 0; i < NEQ; i++)
        {
            // cout << "  ";
            if (i == 0)
                cout << "T";
            else
                cout << "Y" << names[i - 1];
            cout << ": " << J[i][i] << endl;
        }
    }
    if (ToT == 1 || ToT == 2)
    {
        //save to csv files
        ofstream jacFile("sensitivities.csv");
        jacFile << "species";
        jacFile << ",dT_0";
        for (size_t k = 0; k < Nsp; k++)
        {
            jacFile << ",dY" << names[k] << "_0";
        }
        jacFile << endl;
        for (size_t i = 0; i < NEQ; i++)
        {
            if (i == 0)
            {
                jacFile << "dT_f";
            }
            else
            {
                jacFile << "dY" << names[i - 1] << "_f";
            }

            for (size_t j = 0; j < NEQ; j++)
            {
                jacFile << "," << scientific << setprecision(12) << J[i][j];
            }
            jacFile << endl;
        }
        jacFile.close();
    }

    //deallocate
    N_VDestroy(y);
    N_VDestroyVectorArray(yS, Ns);
    CVodeFree(&m_cvode_mem);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
    SUNContext_Free(&sunctx);

    return 0;
}
