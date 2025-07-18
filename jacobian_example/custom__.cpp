#include "cantera/core.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_sens.h>
// #include <sundials/cvodes/cvodes_sens.h>
#include <nvector/nvector_serial.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

namespace fs = std::filesystem;
using namespace Cantera;

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
    /// y: [ T,  Y₁, …, Yₙ ],   ydot same layout
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
    /// fill y0 = [T₀, Y₀₁, …]
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

// CVODE callback: ydot = f(t,y)
static int cvode_rhs(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    ReactorODEs* odes = static_cast<ReactorODEs*>(user_data);
    double* y_data    = NV_DATA_S(y);
    double* yd_data   = NV_DATA_S(ydot);
    odes->eval((double)t, y_data, yd_data);
    return 0;
}

int main() {
    // ——— User settings ———
    std::string mechanism_ = "nDodecane_Reitz.yaml";
    std::string name_      = "nDodecane_IG";
    std::string comp_mix   = "c12h26:1, o2:1, n2:3.76";
    double temp_           = 1001.0;
    double pressure_       = OneAtm;
    const double t0        = 0.0;
    const double tfinal    = 1e-3;
    // dt no longer needed for Jacobian, but kept if you want time‐series later
    const double dt        = 1e-6;

    // — create gas object & set initial state — 
    auto sol = newSolution(mechanism_, name_, "none");
    auto gas = sol->thermo();
    gas->setState_TPX(temp_, pressure_, comp_mix);

    // — build ODE functor —
    ReactorODEs odes(sol);
    size_t NEQ = odes.neq();
    size_t Nsp = odes.nSpecies();
    const auto& names = odes.speciesNames();

    // — allocate CVODE state vector y (size NEQ) & fill y(t0) —
    N_Vector y = N_VNew_Serial(NEQ);
    assert(y);
    double* y_data = NV_DATA_S(y);
    odes.getState(y_data);

    // — create CVODE solver memory and initialize —
    void* cvode_mem = CVodeCreate(CV_BDF);
    assert(cvode_mem);
    int flag = CVodeInit(cvode_mem, cvode_rhs, t0, y);
    assert(flag>=0);
    // scalar tolerances; you can tune these
    double reltol = 1e-8;
    double abstol = 1e-12;
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    assert(flag>=0);
    flag = CVodeSetUserData(cvode_mem, &odes);
    assert(flag>=0);

    // — SENSITIVITY SETUP FOR y0 PARAMETERS — 
    int Ns = (int)NEQ;     // one “parameter” per initial state variable
    N_Vector* yS = N_VCloneVectorArray_Serial(Ns, y);
    assert(yS);
    // S(t0) = I
    for (int is = 0; is < Ns; ++is) {
        for (int j = 0; j < (int)NEQ; ++j) {
            NV_Ith_S(yS[is], j) = (is == j ? 1.0 : 0.0);
        }
    }
    // init forward sensitivities, staggered method
    flag = CVodeSensInit(cvode_mem, Ns, CV_STAGGERED, /*fS*/ nullptr, yS);
    assert(flag>=0);
    flag = CVodeSensEEtolerances(cvode_mem);
    assert(flag>=0);
    // no dependence of f on p beyond initial condition
    flag = CVodeSetSensParams(cvode_mem,
                              /*p*/ nullptr,
                              /*pbar*/ nullptr,
                              /*plist*/ nullptr);
    assert(flag>=0);

    // — INTEGRATE TO tfinal — 
    double t = t0;
    flag = CVode(cvode_mem, tfinal, y, &t, CV_NORMAL);
    assert(flag>=0);

    // — EXTRACT FINAL SENSITIVITIES — 
    flag = CVodeGetSens(cvode_mem, &t, yS);
    assert(flag>=0);

    // — COLLECT INTO A MATRIX J[ i ][ j ] = ∂y_i / ∂y0_j — 
    std::vector<std::vector<double>> J;
    J.resize(NEQ, std::vector<double>(NEQ));
    for (int j = 0; j < (int)NEQ; ++j) {
        double* Sj = NV_DATA_S(yS[j]);
        for (int i = 0; i < (int)NEQ; ++i) {
            J[i][j] = Sj[i];
        }
    }

    // — WRITE SINGLE CSV file bin/J_final.csv — 
    fs::create_directories("bin");
    std::ofstream fout("bin/J_final.csv");
    fout << std::fixed << std::setprecision(8);
    // header
    fout << "state";
    fout << ",T0";
    for (size_t k = 0; k < Nsp; ++k) {
        fout << ",Y0_" << names[k];
    }
    fout << "\n";
    // rows: final‐T and final‐Y_k
    // row 0 = final T
    fout << "T_final";
    for (size_t j = 0; j < NEQ; ++j) {
        fout << "," << J[0][j];
    }
    fout << "\n";
    // rows 1..Nsp = each species
    for (size_t i = 1; i < NEQ; ++i) {
        fout << names[i-1] << "_final";
        for (size_t j = 0; j < NEQ; ++j) {
            fout << "," << J[i][j];
        }
        fout << "\n";
    }
    fout.close();

    // — CLEAN UP — 
    N_VDestroy(y);
    N_VDestroyVectorArray_Serial(yS, Ns);
    CVodeFree(&cvode_mem);

    return 0;
}
