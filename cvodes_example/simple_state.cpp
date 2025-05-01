#include "cantera/core.h"
#include "cantera/numerics/Integrator.h"
#include <iostream>
#include <memory>

using namespace Cantera;

// class to create object to represent the system of odes for the reactor
// from cantera's Funceval to interface with integrator
class ReactorODEs : public FuncEval {
public:
  // constructor initializes thermodynamic and kinetic models from the solution
  // object
  ReactorODEs(std::shared_ptr<Solution> sol)
      : m_gas(sol->thermo()), m_kinetics(sol->kinetics()),
        m_pressure(m_gas->pressure()), m_nSpecies(m_gas->nSpecies()) {
    // allocate space for species mass fractions and temperature
    m_state.resize(m_nSpecies + 1);
  }

  // returns the number of equations (species mass fractions + temperature)
  size_t neq() const override { return m_nSpecies + 1; }

  // evaluates the time derivatives of the system at time t
  void eval(double t, double *y, double *ydot, double *p) override {

    // update the thermodynamic state based on current temperature and mass
    // fractions
    m_gas->setState_TPY(y[m_nSpecies], m_pressure, y);

    // compute net production rates of species
    std::vector<double> omega(m_nSpecies);
    m_kinetics->getNetProductionRates(omega.data());

    // calculate time derivatives of species mass fractions
    for (size_t i = 0; i < m_nSpecies; ++i) {
      ydot[i] = omega[i] * m_gas->molecularWeight(i) / m_gas->density();
    }

    // calculate time derivative of temperature using the energy equation
    double cp = m_gas->cp_mass();
    double hdot = 0.0;
    std::vector<double> h_molar(m_nSpecies);
    m_gas->getPartialMolarEnthalpies(h_molar.data());
    for (size_t i = 0; i < m_nSpecies; ++i) {
      hdot += omega[i] * h_molar[i];
    }
    ydot[m_nSpecies] = -hdot / (m_gas->density() * cp);
  }

  // provides the ic (initial species mass fractions and intial temp.) for the
  // integrator
  void getInitialConditions(double t0, double *y) {
    m_gas->getMassFractions(y);
    y[m_nSpecies] = m_gas->temperature();
  }

  // added implementation for getstate() to avoid notimplementederror.
  void getState(double *y) override {
    m_gas->getMassFractions(y);
    y[m_nSpecies] = m_gas->temperature();
  }

private:
  std::shared_ptr<ThermoPhase> m_gas;
  std::shared_ptr<Kinetics> m_kinetics;
  double m_pressure;
  size_t m_nSpecies;
  std::vector<double>
      m_state; // state vector: species mass fractions and temperature
};

int main() {
  // create an object
  auto sol = newSolution("gri30.yaml");

  // access the thermodynamic model and set the initial state
  auto gas = sol->thermo();
  gas->setState_TPX(1000.0, OneAtm, "H2:2,O2:1,N2:4");

  // create new instace of reactorodes object with the solution
  ReactorODEs reactor(sol);

  // create and initialize the cvode integrator
  auto integrator = newIntegrator("CVODE");
  integrator->initialize(0.0, reactor); // 
  integrator->setTolerances(1e-9, 1e-15);

  // integrate the system up to the specified end time
  double t_end = 1.0;
  integrator->integrate(t_end);

  // retrieve the solution vector after integration
  double *y = integrator->solution();

  // print results
  std::cout << y[gas->nSpecies()] << " [K]" << std::endl;
  delete integrator;

  return 0;
}
