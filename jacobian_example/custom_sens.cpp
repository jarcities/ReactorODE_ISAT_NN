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

//------------------------------------------------------------------------------
class ReactorODEs {
public:
    ReactorODEs(shared_ptr<Solution> sol) {
        m_gas      = sol->thermo();
        m_kinetics = sol->kinetics();
        m_pressure = m_gas->pressure();
        m_nSpecies = m_gas->nSpecies();
        m_nEqs     = m_nSpecies + 1;
        m_hbar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);
    }

    void eval(double /*t*/, const double* y, double* ydot) {
        double T = y[0];
        const double* Y = y + 1;
        m_gas->setMassFractions_NoNorm(Y);
        m_gas->setState_TP(T, m_pressure);
        double rho = m_gas->density();
        double cp  = m_gas->cp_mass();
        m_gas->getPartialMolarEnthalpies(m_hbar.data());
        m_kinetics->getNetProductionRates(m_wdot.data());
        
        // dT/dt
        double hdot_vol = 0.0;
        for (size_t k = 0; k < m_nSpecies; k++) {
            hdot_vol += m_hbar[k] * m_wdot[k];
        }
        ydot[0] = - hdot_vol / (rho * cp);
        
        // dYk/dt
        for (size_t k = 0; k < m_nSpecies; k++) {
            ydot[1 + k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
        }
    }

    void getState(double* y) {
        y[0] = m_gas->temperature();
        m_gas->getMassFractions(y + 1);
    }
    
    size_t neq() const {
        return m_nEqs;
    }
    
    size_t nSpecies() const {
        return m_nSpecies;
    }
    
    const std::vector<std::string>& speciesNames() const {
        return m_gas->speciesNames();
    }

private:
    shared_ptr<ThermoPhase> m_gas;
    shared_ptr<Kinetics>    m_kinetics;
    double                  m_pressure;
    size_t                  m_nSpecies, m_nEqs;
    std::vector<double>     m_hbar, m_wdot;
};

static int cvodes_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
    ReactorODEs* odes = static_cast<ReactorODEs*>(user_data);
    double* y_data    = NV_DATA_S(y);
    double* yd_data   = NV_DATA_S(ydot);
    odes->eval((double)t, y_data, yd_data);
    return 0;
}

int main() {
    std::string mechanism_ = "h2o2.yaml";
    std::string name_      = "ohmech";
    std::string comp_mix   = "H2:2, O2:1, N2:4";
    double temp_           = 1001.0;
    double pressure_       = OneAtm;
    const double t0        = 0.0;
    const double tfinal    = 1e-4;  

    //init gas object and IC
    auto sol = newSolution(mechanism_, name_, "none");
    auto gas = sol->thermo();
    gas->setState_TPX(temp_, pressure_, comp_mix);

    //init reactor object
    ReactorODEs odes(sol);
    size_t NEQ = odes.neq();
    size_t Nsp = odes.nSpecies();
    const auto& names = odes.speciesNames();

    cout << "\nInitial conditions:" << endl;
    cout << "Temperature: " << temp_ << " K" << endl;
    for (size_t k = 0; k < Nsp; k++) {
        double Yk = gas->massFraction(k);  // Use massFraction instead of getMassFraction
        cout << "Y_" << names[k] << ": " << Yk << endl;
    }
    cout << "Number of equations: " << NEQ << endl;

    // Create SUNDIALS context (required for SUNDIALS 6.0+)
    SUNContext sunctx;
    int flag = SUNContext_Create(SUN_COMM_NULL, &sunctx);
    assert(flag == 0);

    // — allocate CVODES state vector y (size NEQ) & fill y(t0) —
    N_Vector y = N_VNew_Serial(NEQ, sunctx);
    assert(y);
    double* y_data = NV_DATA_S(y);
    odes.getState(y_data);

    // — create CVODES solver memory and initialize —
    void* cvodes_mem = CVodeCreate(CV_BDF, sunctx);
    assert(cvodes_mem);
    flag = CVodeInit(cvodes_mem, cvodes_rhs, t0, y);
    assert(flag >= 0);
    
    // scalar tolerances - relaxed for stiff chemical kinetics
    double reltol = 1e-6;
    double abstol = 1e-9;
    flag = CVodeSStolerances(cvodes_mem, reltol, abstol);
    assert(flag >= 0);
    flag = CVodeSetUserData(cvodes_mem, &odes);
    assert(flag >= 0);

    // Set max number of steps to handle stiff problems
    flag = CVodeSetMaxNumSteps(cvodes_mem, 50000);
    assert(flag >= 0);

    // Set max step size for stability
    flag = CVodeSetMaxStep(cvodes_mem, 1e-6);
    assert(flag >= 0);

    // — linear solver —
    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    assert(A);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    assert(LS);
    flag = CVodeSetLinearSolver(cvodes_mem, LS, A);
    assert(flag >= 0);

    cout << "\nSetting up sensitivity analysis..." << endl;

    // — SENSITIVITY SETUP FOR y0 PARAMETERS —
    int Ns = (int)NEQ;     // one "parameter" per initial state variable
    N_Vector* yS = N_VCloneVectorArray(Ns, y);
    assert(yS);
    
    // Initialize sensitivity matrix S(t0) = I (identity matrix)
    // This means we're computing dy/dy0 where y0 are the initial conditions
    for (int is = 0; is < Ns; ++is) {
        for (int j = 0; j < (int)NEQ; ++j) {
            NV_Ith_S(yS[is], j) = (is == j ? 1.0 : 0.0);
        }
    }
    
    // Initialize forward sensitivities using staggered method
    // NULL for fS means CVODES will compute sensitivities automatically
    flag = CVodeSensInit(cvodes_mem, Ns, CV_STAGGERED, /*fS*/ nullptr, yS);
    assert(flag >= 0);
    
    // Use estimated error tolerances for sensitivities
    flag = CVodeSensEEtolerances(cvodes_mem);
    assert(flag >= 0);

    // Set up parameter array for sensitivity analysis
    // We need to provide parameter values even though we're computing dy/dy0
    vector<sunrealtype> p(NEQ);
    vector<sunrealtype> pbar(NEQ);
    for (int i = 0; i < NEQ; i++) {
        p[i] = y_data[i];      // Initial conditions as parameters
        pbar[i] = (fabs(y_data[i]) > 1e-12) ? fabs(y_data[i]) : 1.0;  // Scaling factors
    }

    // Set sensitivity parameters
    flag = CVodeSetSensParams(cvodes_mem,
                              p.data(),        // parameter values
                              pbar.data(),     // parameter scaling factors  
                              nullptr);        // parameter list (NULL means all)
    assert(flag >= 0);

    cout << "Number of sensitivity parameters: " << Ns << endl;
    cout << "Sensitivity method: CV_STAGGERED" << endl;

    // — INTEGRATE TO tfinal —
    cout << "\nIntegrating from t = " << t0 << " to t = " << tfinal << " s..." << endl;
    double t = t0;
    flag = CVode(cvodes_mem, tfinal, y, &t, CV_NORMAL);
    
    if (flag < 0) {
        cout << "CVode integration failed with flag = " << flag << endl;
        return 1;
    }
    
    cout << "Integration successful! Final time: " << t << " s" << endl;

    // Get final solution
    cout << "\nFinal conditions:" << endl;
    cout << "Temperature: " << y_data[0] << " K" << endl;
    for (size_t k = 0; k < Nsp; k++) {
        cout << "Y_" << names[k] << ": " << y_data[k+1] << endl;
    }

    // — EXTRACT FINAL SENSITIVITIES —
    cout << "\nExtracting final sensitivity matrix..." << endl;
    flag = CVodeGetSens(cvodes_mem, &t, yS);
    assert(flag >= 0);

    // — COLLECT INTO A MATRIX J[i][j] = ∂y_final[i] / ∂y_initial[j] —
    std::vector<std::vector<double>> J;
    J.resize(NEQ, std::vector<double>(NEQ));
    for (int j = 0; j < (int)NEQ; ++j) {
        double* Sj = NV_DATA_S(yS[j]);
        for (int i = 0; i < (int)NEQ; ++i) {
            J[i][j] = Sj[i];
        }
    }

    // — PRINT JACOBIAN MATRIX —
    cout << "\nJacobian Matrix dy_final/dy_initial:" << endl;
    cout << "====================================" << endl;
    
    // Print header
    cout << setw(12) << " ";
    cout << setw(15) << "dT_0";
    for (size_t k = 0; k < Nsp; k++) {
        cout << setw(15) << ("dY" + names[k] + "_0");
    }
    cout << endl;
    
    // Print rows
    for (size_t i = 0; i < NEQ; i++) {
        if (i == 0) {
            cout << setw(12) << "dT_f:";
        } else {
            cout << setw(12) << ("dY" + names[i-1] + "_f:");
        }
        
        for (size_t j = 0; j < NEQ; j++) {
            cout << setw(15) << scientific << setprecision(6) << J[i][j];
        }
        cout << endl;
    }

    // — SAVE JACOBIAN TO FILE —
    ofstream jacFile("jacobian.csv");
    
    //header
    jacFile << "Variable";
    jacFile << ",dT_0";
    for (size_t k = 0; k < Nsp; k++) {
        jacFile << ",dY" << names[k] << "_0";
    }
    jacFile << endl;
    
    //jacobian
    for (size_t i = 0; i < NEQ; i++) {
        if (i == 0) {
            jacFile << "dT_f";
        } else {
            jacFile << "dY" << names[i-1] << "_f";
        }
        
        for (size_t j = 0; j < NEQ; j++) {
            jacFile << "," << scientific << setprecision(12) << J[i][j];
        }
        jacFile << endl;
    }
    jacFile.close();

    // — ANALYZE JACOBIAN —
    cout << "\nJacobian Analysis:" << endl;
    cout << "==================" << endl;
    
    // Find largest elements
    double max_element = 0.0;
    size_t max_i = 0, max_j = 0;
    
    for (size_t i = 0; i < NEQ; i++) {
        for (size_t j = 0; j < NEQ; j++) {
            if (fabs(J[i][j]) > fabs(max_element)) {
                max_element = J[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    cout << "Largest Jacobian element: " << max_element << endl;
    cout << "Location: ";
    if (max_i == 0) cout << "dT_f"; 
    else cout << "dY" << names[max_i-1] << "_f";
    cout << " / ";
    if (max_j == 0) cout << "dT_0"; 
    else cout << "dY" << names[max_j-1] << "_0";
    cout << endl;
    
    // Diagonal elements (sensitivity of each variable to itself)
    cout << "\nDiagonal elements (self-sensitivities):" << endl;
    for (size_t i = 0; i < NEQ; i++) {
        cout << "  ";
        if (i == 0) cout << "T";
        else cout << "Y" << names[i-1];
        cout << ": " << J[i][i] << endl;
    }

    //deallocate
    N_VDestroy(y);
    N_VDestroyVectorArray(yS, Ns);
    CVodeFree(&cvodes_mem);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
    SUNContext_Free(&sunctx);

    return 0;
}
