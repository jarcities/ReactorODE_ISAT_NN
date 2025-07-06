////////////////////////
////ORIGINAL EXAMPLE////
////////////////////////
// /*
//  * Custom reactor
//  * ==============
//  *
//  * Solve a closed-system constant pressure ignition problem where the governing
//  * equations are custom-implemented, using Cantera's interface to CVODES to
//  * integrate the equations.
//  *
//  * .. tags:: C++, combustion, reactor network, user-defined model, ignition delay
//  */

// // This file is part of Cantera. See License.txt in the top-level directory or
// // at https://cantera.org/license.txt for license and copyright information.

// #include "cantera/core.h"
// #include "cantera/numerics/Integrator.h"
// #include <fstream>
// #include <iostream>

// using namespace Cantera;

// class ReactorODEs : public FuncEval {
// public:
//     /**
//      * Constructor
//      * @param[in] sol Solution object specifying initial system state.
//      */
//     ReactorODEs(shared_ptr<Solution> sol) {
//         /* ---------------------- INITIALIZE MEMBER VARS ---------------------- */

//         // pointer to the system's ThermoPhase object. updated by the solver during
//         // simulation to provide iteration-specific thermodynamic properties.
//         m_gas = sol->thermo();

//         // pointer to the kinetics manager. provides iteration-specific species
//         // production rates based on the current state of the ThermoPhase.
//         m_kinetics = sol->kinetics();

//         // the system's constant pressure, taken from the provided initial state.
//         m_pressure = m_gas->pressure();

//         // number of chemical species in the system.
//         m_nSpecies = m_gas->nSpecies();

//         // resize the vector<double> storage containers for species partial molar enthalpies
//         // and net production rates. internal values are updated and used by the solver
//         // per iteration.
//         m_hbar.resize(m_nSpecies);
//         m_wdot.resize(m_nSpecies);

//         // number of equations in the ODE system. a conservation equation for each
//         // species, plus a single energy conservation equation for the system.
//         m_nEqs = m_nSpecies + 1;
//     }

//     /* %%
//      * Evaluate the ODE right-hand-side function, :math:`\dot{y} = f(t,y)`.
//      *
//      * Overridden from :ct:`FuncEval`, called by the integrator during simulation.
//      *
//      * :param t:
//      *     time
//      * :param y:
//      *     solution vector, length neq()
//      * :param ydot:
//      *     rate of change of solution vector, length neq()
//      * :param p:
//      *     sensitivity parameter vector, length nparams()
//      *
//      * note: sensitivity analysis isn't implemented in this example
//      */
//     void eval(double t, double* y, double* ydot, double* p) override {
//         // the solution vector *y* is [T, Y_1, Y_2, ... Y_K], where T is the
//         // system temperature, and Y_k is the mass fraction of species k.
//         // similarly, the time derivative of the solution vector, *ydot*, is
//         // [dT/dt, Y_1/dt, Y_2/dt, ... Y_K/dt].
//         // the following variables are defined for clear and convenient access
//         // to these vectors:
//         double temperature = y[0];
//         double *massFracs = &y[1];
//         double *dTdt = &ydot[0];
//         double *dYdt = &ydot[1];

//         /* ------------------------- UPDATE GAS STATE ------------------------- */
//         // the state of the ThermoPhase is updated to reflect the current solution
//         // vector, which was calculated by the integrator.
//         m_gas->setMassFractions_NoNorm(massFracs);
//         m_gas->setState_TP(temperature, m_pressure);

//         /* ----------------------- GET REQ'D PROPERTIES ----------------------- */
//         double rho = m_gas->density();
//         double cp = m_gas->cp_mass();
//         m_gas->getPartialMolarEnthalpies(&m_hbar[0]);
//         m_kinetics->getNetProductionRates(&m_wdot[0]);

//         /* -------------------------- ENERGY EQUATION ------------------------- */
//         // the rate of change of the system temperature is found using the energy
//         // equation for a closed-system constant pressure ideal gas:
//         //     m*cp*dT/dt = - sum[h(k) * dm(k)/dt]
//         // or equivalently:
//         //     dT/dt = - sum[hbar(k) * dw(k)/dt] / (rho * cp)
//         double hdot_vol = 0;
//         for (size_t k = 0; k < m_nSpecies; k++) {
//             hdot_vol += m_hbar[k] * m_wdot[k];
//         }
//         *dTdt = - hdot_vol / (rho * cp);

//         /* --------------------- SPECIES CONSERVATION EQS --------------------- */
//         // the rate of change of each species' mass fraction is found using the closed-system
//         // species conservation equation, applied once for each species:
//         //     m*dY(k)/dt = dm(k)/dt
//         // or equivalently:
//         //     dY(k)/dt = dw(k)/dt * MW(k) / rho
//         for (size_t k = 0; k < m_nSpecies; k++) {
//             dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
//         }
//     }

//     /**
//      * Number of equations in the ODE system.
//      *   - overridden from FuncEval, called by the integrator during initialization.
//      */
//     size_t neq() const override {
//         return m_nEqs;
//     }

//     /**
//      * Provide the current values of the state vector, *y*.
//      *   - overridden from FuncEval, called by the integrator during initialization.
//      * @param[out] y solution vector, length neq()
//      */
//     void getState(double* y) override {
//         // the solution vector *y* is [T, Y_1, Y_2, ... Y_K], where T is the
//         // system temperature, and Y_k is the mass fraction of species k.
//         y[0] = m_gas->temperature();
//         m_gas->getMassFractions(&y[1]);
//     }

// private:
//     // private member variables, to be used internally.
//     shared_ptr<ThermoPhase> m_gas;
//     shared_ptr<Kinetics> m_kinetics;
//     vector<double> m_hbar;
//     vector<double> m_wdot;
//     double m_pressure;
//     size_t m_nSpecies;
//     size_t m_nEqs;
// };

// int main() {
//     /* -------------------- CREATE GAS & SPECIFY STATE -------------------- */
//     auto sol = newSolution("h2o2.yaml", "ohmech", "none");
//     auto gas = sol->thermo();
//     gas->setState_TPX(1001, OneAtm, "H2:2, O2:1, N2:4");

//     /* ---------------------- CREATE CSV OUTPUT FILE ---------------------- */
//     // simulation results will be outputted to a .csv file as complete state vectors
//     // for each simulation time point.
//     // create the csv file, overwriting any existing ones with the same name.
//     std::std::ofstream outputFile("custom_cxx.csv");

//     // for convenience, a title for each of the state vector's components is written to
//     // the first line of the csv file.
//     outputFile << "time (s), temp (K)";
//     for (size_t k = 0; k < gas->nSpecies(); k++) {
//         outputFile << ", Y_" << gas->speciesName(k);
//     }
//     outputFile << std::endl;

//     /* --------------------- CREATE ODE RHS EVALUATOR --------------------- */
//     ReactorODEs odes = ReactorODEs(sol);

//     /* ---------------------- SPECIFY TIME INTERVAL ----------------------- */
//     // the simulation is run over the time interval specified below. tnow is initialized
//     // with the simulation's start time, and is updated on each timestep to reflect
//     // the new absolute time of the system.
//     double tnow = 0.0;
//     double tfinal = 1e-3;

//     /* ------------------- CREATE & INIT ODE INTEGRATOR ------------------- */
//     // create an ODE integrator object, which will be used to solve the system of ODES defined
//     // in the ReactorODEs class. a C++ interface to the C-implemented SUNDIALS CVODES integrator
//     // (CVodesIntegrator) is built into Cantera, and will be used to solve this example.
//     //  - the default settings for CVodesIntegrator are used:
//     //     solution method: BDF_Method
//     //     problem type: DENSE + NOJAC
//     //     relative tolerance: 1.0e-9
//     //     absolute tolerance: 1.0e-15
//     //     max step size: +inf
//     unique_ptr<Integrator> integrator(newIntegrator("CVODE"));

//     // initialize the integrator, specifying the start time and the RHS evaluator object.
//     // internally, the integrator will apply settings, allocate needed memory, and populate
//     // this memory with the appropriate initial values for the system.
//     integrator->initialize(tnow, odes);

//     /* ----------------------- SIMULATION TIME LOOP ----------------------- */
//     while (tnow < tfinal) {
//         // advance the simulation to the current absolute time, tnow, using the integrator's
//         // ODE system time-integration capability. a pointer to the resulting system state vector
//         // is fetched in order to access the solution components.
//         integrator->integrate(tnow);
//         double *solution = integrator->solution();

//         // output the current absolute time and solution state vector to the csv file. you can view
//         // these results by opening the "custom_cxx.csv" file that appears in this program's directory
//         // after compiling and running.
//         outputFile << tnow;
//         for (size_t i = 0; i < odes.neq(); i++) {
//             outputFile << ", " << solution[i];
//         }
//         outputFile << std::endl;

//         // increment the simulation's absolute time, tnow, then return to the start of the loop to
//         // advance the simulation to this new time point.
//         tnow += 1e-5;
//     }
// }













































// ///////////////////////////
// ////USING SENSITIVITIES////
// ///////////////////////////
// #include "cantera/core.h"
// #include "cantera/numerics/Integrator.h"
// #include <fstream>
// #include <iostream>
// #include <filesystem>
// #include <algorithm>
// namespace fs = std::filesystem;

// using namespace Cantera;

// class ReactorODEs : public FuncEval
// {
// public:
//     ReactorODEs(shared_ptr<Solution> sol)
//     {
//         m_gas = sol->thermo();
//         m_kin = sol->kinetics();
//         m_pressure = m_gas->pressure();
//         m_nSpecies = m_gas->nSpecies();

//         // allocate storage for enthalpies and reaction rates
//         m_hbar.resize(m_nSpecies);
//         m_wdot.resize(m_nSpecies);
//         // m_nEqs = m_nSpecies + 1; //this
//         m_nPrimary = m_nSpecies + 1; // that
//         m_nEqs = m_nPrimary * 2;     // that

//         // sensitivity calc, only using one parameter (temp.)
//         m_sens_params.resize(1);                 // sensitivity calc
//         m_paramScales.resize(1);                 // sensitivity calc
//         m_sens_params[0] = m_gas->temperature(); // sensitivity calc
//         m_paramScales[0] = m_gas->temperature(); // sensitivity calc
//     }

//     // void eval(double t, double *y, double *ydot, double *p) override
//     // {
//     //     double T = y[0];
//     //     double *Ys = &y[1];
//     //     double *dTdt = &ydot[0];
//     //     double *dYdt = &ydot[1];

//     //     m_gas->setMassFractions_NoNorm(Ys);
//     //     m_gas->setState_TP(T, m_pressure);

//     //     double rho = m_gas->density();
//     //     double cp = m_gas->cp_mass();
//     //     m_gas->getPartialMolarEnthalpies(m_hbar.data());
//     //     m_kin->getNetProductionRates(m_wdot.data());

//     //     double hdot_vol = 0.0;
//     //     for (size_t k = 0; k < m_nSpecies; ++k)
//     //     {
//     //         hdot_vol += m_hbar[k] * m_wdot[k];
//     //     }
//     //     *dTdt = -hdot_vol / (rho * cp);

//     //     for (size_t k = 0; k < m_nSpecies; ++k)
//     //     {
//     //         dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
//     //     }

//     //     //............................................................
//     //     double sT = y[m_nPrimary + 0];

//     //     vector<double> dwdot_dT(m_nSpecies);
//     //     m_kin->getNetProductionRates_ddT(dwdot_dT.data());

//     //     double eps = 1e-8 * std::max(1.0, T);
//     //     vector<double> hbar_p(m_nSpecies), hbar_m(m_nSpecies);
//     //     m_gas->setState_TP(T + eps, m_pressure);
//     //     m_gas->getPartialMolarEnthalpies(hbar_p.data());
//     //     m_gas->setState_TP(T - eps, m_pressure);
//     //     m_gas->getPartialMolarEnthalpies(hbar_m.data());
//     //     m_gas->setState_TP(T, m_pressure);

//     //     double sum1 = 0.0, sum2 = 0.0;
//     //     for (size_t k = 0; k < m_nSpecies; ++k)
//     //     {
//     //         double dhbar_dT = (hbar_p[k] - hbar_m[k]) / (2 * eps);
//     //         sum1 += m_hbar[k] * dwdot_dT[k];
//     //         sum2 += m_wdot[k] * dhbar_dT;
//     //     }
//     //     double df_dT = -(sum1 + sum2) / (rho * cp);
//     //     ydot[m_nPrimary + 0] = df_dT * sT;

//     //     double drho_dT = -rho * m_gas->thermalExpansionCoeff();
//     //     for (size_t k = 0; k < m_nSpecies; ++k)
//     //     {
//     //         double MW = m_gas->molecularWeight(k);
//     //         double term = dwdot_dT[k] + m_wdot[k] * (drho_dT / rho);
//     //         ydot[m_nPrimary + 1 + k] = (MW / rho) * term * sT;
//     //     }
//     //     //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     // }
//     void eval(double t, double* y, double* ydot, double* p) override
//     {
//         double T = y[0];
//         double* Ys  = &y[1];
//         double* dTdt = &ydot[0];
//         double* dYdt = &ydot[1];

//         m_gas->setMassFractions_NoNorm(Ys);
//         m_gas->setState_TP(T, m_pressure);

//         double rho = m_gas->density();
//         double cp  = m_gas->cp_mass();
//         m_gas->getPartialMolarEnthalpies(m_hbar.data());
//         m_kin->getNetProductionRates       (m_wdot.data());

//         double hdot_vol = 0.0;
//         for (size_t k = 0; k < m_nSpecies; ++k) {
//             hdot_vol += m_hbar[k] * m_wdot[k];
//         }
//         *dTdt = -hdot_vol / (rho * cp);

//         for (size_t k = 0; k < m_nSpecies; ++k) {
//             double MW = m_gas->molecularWeight(k);
//             dYdt[k] = m_wdot[k] * MW / rho;
//         }

//         //............................................................
//         double sT = y[m_nPrimary + 0];

//         vector<double> dwdot_dT(m_nSpecies);
//         m_kin->getNetProductionRates_ddT(dwdot_dT.data());

//         double eps = 1e-5 * std::max(1.0, T);
//         vector<double> hbar_p(m_nSpecies), hbar_m(m_nSpecies);
//         m_gas->setState_TP(T + eps, m_pressure);
//         m_gas->getPartialMolarEnthalpies(hbar_p.data());
//         m_gas->setState_TP(T - eps, m_pressure);
//         m_gas->getPartialMolarEnthalpies(hbar_m.data());
//         m_gas->setState_TP(T, m_pressure);

//         m_gas->setState_TP(T + eps, m_pressure);
//         double cp_p = m_gas->cp_mass();
//         m_gas->setState_TP(T - eps, m_pressure);
//         double cp_m = m_gas->cp_mass();
//         m_gas->setState_TP(T, m_pressure);
//         double d_cp_dT = (cp_p - cp_m) / (2 * eps);

//         double d_rho_dT = -rho * m_gas->thermalExpansionCoeff();

//         double sum1 = 0.0, sum2 = 0.0;
//         for (size_t k = 0; k < m_nSpecies; ++k) {
//             double dh_dT = (hbar_p[k] - hbar_m[k]) / (2 * eps);
//             sum1 += m_hbar[k] * dwdot_dT[k];
//             sum2 += m_wdot[k]  * dh_dT;
//         }
//         double df_dT = -(sum1 + sum2)/(rho*cp)
//              + hdot_vol * ((rho*d_cp_dT + cp*d_rho_dT)
//                            /((rho*cp)*(rho*cp)));
//         ydot[m_nPrimary + 0] = df_dT * sT;

//         for (size_t k = 0; k < m_nSpecies; ++k) {
//             double MW     = m_gas->molecularWeight(k);
//             double term   = dwdot_dT[k] + m_wdot[k] * (d_rho_dT / rho);
//             ydot[m_nPrimary + 1 + k] = (MW / rho) * term * sT;
//         }
//         //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//     }

//     size_t neq() const override
//     {
//         return m_nEqs;
//     }

//     void getState(double *y) override
//     {
//         y[0] = m_gas->temperature();
//         m_gas->getMassFractions(&y[1]);

//         // sensitivity calc
//         y[m_nPrimary + 0] = 1.0;
//         for (size_t k = 1; k < m_nPrimary; ++k)
//         {
//             y[m_nPrimary + k] = 0.0;
//         }
//     }

//     // sensitivity calc
//     size_t nparams() const override
//     {
//         return 1;
//     }

//     // sensitivity calc
//     void getSensParams(double *p) const
//     {
//         p[0] = m_sens_params[0];
//     }

//     // sensitivity calc
//     void getSensScales(double *p) const
//     {
//         p[0] = m_paramScales[0];
//     }

// private:
//     shared_ptr<ThermoPhase> m_gas;
//     shared_ptr<Kinetics> m_kin;
//     vector<double> m_hbar;
//     vector<double> m_wdot;
//     double m_pressure;
//     size_t m_nSpecies;
//     size_t m_nEqs;
//     size_t m_nPrimary;
// };

// int main()
// {
//     // init new gas
//     auto sol = newSolution("h2o2.yaml", "ohmech", "none");
//     auto gas = sol->thermo();
//     gas->setState_TPX(1001, OneAtm, "H2:2, O2:1, N2:4");

//     // new reactor object
//     ReactorODEs odes(sol);
//     double t0 = 0.0;
//     double tfinal = 1e-3;
//     double dt = 1e-5;
//     double reltol = 1e-9;
//     double abstol = 1e-15;
//     double sensreltol = 1e-9;
//     double sensabstol = 1e-15;

//     // cvodes object and functions
//     auto integrator = newIntegrator("CVODE");
//     integrator->setMethod(BDF_Method);        // sensitivity calc
//     integrator->setLinearSolverType("DENSE"); // sensitivity calc
//     integrator->setTolerances(reltol, abstol);
//     integrator->setSensitivityTolerances(sensreltol, sensabstol); // sensitivity calc
//     integrator->initialize(t0, odes);
//     size_t npar = odes.nparams();
//     std::vector<double> s0(odes.neq(), 0.0);
//     s0[0] = 1.0;

//     //..........................................................
//     // dt for central difference
//     const double dT0 = 1e-5;

//     // forward perturbed
//     auto sol_f = newSolution("h2o2.yaml", "ohmech", "none");
//     auto gas_f = sol_f->thermo();
//     gas_f->setState_TPX(1001 + dT0, OneAtm, "H2:2, O2:1, N2:4");
//     ReactorODEs odes_f(sol_f);
//     auto int_f = newIntegrator("CVODE");
//     int_f->setMethod(BDF_Method);
//     int_f->setLinearSolverType("DENSE");
//     int_f->setTolerances(reltol, abstol);
//     int_f->initialize(t0, odes_f);

//     // backward perturbed
//     auto sol_b = newSolution("h2o2.yaml", "ohmech", "none");
//     auto gas_b = sol_b->thermo();
//     gas_b->setState_TPX(1001 - dT0, OneAtm, "H2:2, O2:1, N2:4");
//     ReactorODEs odes_b(sol_b);
//     auto int_b = newIntegrator("CVODE");
//     int_b->setTolerances(reltol, abstol);
//     int_b->setMethod(BDF_Method);
//     int_b->setLinearSolverType("DENSE");
//     int_b->setTolerances(reltol, abstol);
//     int_b->initialize(t0, odes_b);
//     //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//     // create csv file headers
//     fs::create_directories("bin");
//     std::vector<std::std::ofstream> file(gas->nSpecies());
//     for (size_t k = 0; k < gas->nSpecies(); ++k)
//     {
//         std::string fname = "bin/" + gas->speciesName(k) + ".csv";
//         file[k].open(fname);
//         file[k] << "time,temp,"
//                 << "Y_" << gas->speciesName(k) << ','
//                 << "dY_" << gas->speciesName(k) << "/dT0,"
//                 << "dY_" << gas->speciesName(k) << "/dT0_CD\n";
//     }

//     double t = t0;
//     while (t < tfinal)
//     {
//         t += dt;
//         integrator->integrate(t);
//         double *soln = integrator->solution();

//         //..............................
//         // central difference
//         int_f->integrate(t);
//         int_b->integrate(t);
//         double *yf = int_f->solution();
//         double *yb = int_b->solution();
//         //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//         for (size_t k = 0; k < gas->nSpecies(); ++k)
//         {
//             double *soln = integrator->solution();
//             double *sens = soln + odes.nparams() * 0 + odes.neq() / 2;
//             file[k] << t
//                     << ',' << soln[0]
//                     << ',' << soln[1 + k]
//                     << ',' << sens[1 + k]
//                     << ',' << (yf[1 + k] - yb[1 + k]) / (2.0 * dT0)
//                     << '\n';
//         }
//     }

//     return 0;
// }














































/////////////////////
////WITH JACOBIAN////
/////////////////////
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

namespace fs = std::filesystem;
using namespace Cantera;

int main()
{

    //timing
    std::chrono::duration<double> jac_total{};
    std::chrono::duration<double> cd_total{};

    //init solution and states
    auto sol = newSolution("h2o2.yaml", "ohmech", "none");
    auto gas = sol->thermo();
    gas->setState_TPX(1001.0, OneAtm, "H2:2, O2:1, N2:4");
    const double t0 = 0.0;
    const double tfinal = 1e-3;
    const double dt = 1e-5;
    const double dT0 = 1e-5;

    //new reactor object
    auto base_r = newReactor("ConstPressureReactor", sol);
    auto reactor = std::static_pointer_cast<Reactor>(base_r);
    ReactorNet net;
    net.addReactor(*reactor);

    //forward perturbation for central difference
    auto sol_f = newSolution("h2o2.yaml", "ohmech", "none");
    auto gas_f = sol_f->thermo();
    gas_f->setState_TPX(1001.0 + dT0, OneAtm, "H2:2, O2:1, N2:4");
    auto base_rf = newReactor("ConstPressureReactor", sol_f);
    auto r_f = std::static_pointer_cast<Reactor>(base_rf);
    ReactorNet net_f;
    net_f.addReactor(*r_f);

    //backward perturbation for central difference
    auto sol_b = newSolution("h2o2.yaml", "ohmech", "none");
    auto gas_b = sol_b->thermo();
    gas_b->setState_TPX(1001.0 - dT0, OneAtm, "H2:2, O2:1, N2:4");
    auto base_rb = newReactor("ConstPressureReactor", sol_b);
    auto r_b = std::static_pointer_cast<Reactor>(base_rb);
    ReactorNet net_b;
    net_b.addReactor(*r_b);

    //init all ReactorNet objects
    net.setInitialTime(t0);
    net.initialize();
    net_f.setInitialTime(t0);
    net_f.initialize();
    net_b.setInitialTime(t0);
    net_b.initialize();

    //open csv files
    size_t n_species = gas->nSpecies();
    fs::create_directories("bin");
    std::vector<std::ofstream> files(n_species);
    for (size_t k = 0; k < n_species; k++)
    {
        auto name = gas->speciesName(k);
        files[k].open("bin/" + name + ".csv");
        files[k]
            << "time,temp,"
            << "Y_" << name << ','
            << "finiteDifferenceJacobian(),"
            << "reg. central difference\n";
        files[k] << std::fixed << std::setprecision(8);
    }

    //init state vectors
    size_t neq = reactor->neq();
    std::vector<double> y(neq), yf(neq), yb(neq);

    //INTEGRATION LOOP
    for (double t = t0 + dt; t <= tfinal + 1e-12; t += dt)
    {
        //integrate solution
        net.advance(t);
        net_f.advance(t); //cent_diff
        net_b.advance(t); //cent_diff

        //get state solutions
        reactor->getState(y.data());
        r_f->getState(yf.data());
        r_b->getState(yb.data());

        //calc eigen jacobian
        auto jac_start = std::chrono::high_resolution_clock::now(); //timing
        Eigen::SparseMatrix<double> eig_jac_sparse = reactor->finiteDifferenceJacobian();
        Eigen::MatrixXd eig_jac = Eigen::MatrixXd(eig_jac_sparse); //convert back to dense
        auto jac_end = std::chrono::high_resolution_clock::now(); //timing
        std::chrono::duration<double> jac_diff = jac_end - jac_start; //timting
        jac_total += jac_diff; //timing

        //go through each species
        for (size_t k = 0; k < n_species; k++)
        {

            //state vectors
            std::vector<double> y0(neq), ytmp(neq), ydot_f(neq), ydot_b(neq);

            //get state
            reactor->getState(y0.data());
            ytmp = y0;

            //calc central difference at each species
            auto cd_start = std::chrono::high_resolution_clock::now(); //timing
            ytmp[0] = y0[0] + dT0;
            net.eval(t, ytmp.data(), ydot_f.data(), nullptr);
            ytmp[0] = y0[0] - dT0;
            net.eval(t, ytmp.data(), ydot_b.data(), nullptr);
            double central_diff = (ydot_f[1 + k] - ydot_b[1 + k]) / (2.0 * dT0);
            double Yk = y[1 + k];
            auto cd_end = std::chrono::high_resolution_clock::now(); //timing
            std::chrono::duration<double> cd_diff = cd_end - cd_start; //timing
            cd_total += cd_diff; //timing

            //get eigen jacobian at species
            double jac_fd = eig_jac(k + 1, 0);

            //print solutions
            files[k]
                << t << ','
                << y[0] << ','
                << Yk << ','
                << jac_fd << ','
                << central_diff << '\n';
        }
    }

    //print timing
    std::cout << "eigen jacobian time = " << jac_total.count() << " seconds" << std::endl;
    std::cout << "central difference time = " << cd_total.count() << " seconds" << std::endl;

    return 0;
}
