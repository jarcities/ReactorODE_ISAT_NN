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
#include <string>
#include <iomanip>

namespace fs = std::filesystem;
using namespace Cantera;

int main()
{

    //timing
    std::chrono::duration<double> jac_total{};
    std::chrono::duration<double> cd_total{};

    //set mechanisms
    std::string mechanism_ = "nDodecane_Reitz.yaml";
    std::string name_ = "nDodecane_IG";
    std::string comp_mix = "c12h26:1, o2:1, n2:3.76";
    double temp_ = 1001.0; //kelvin
    double pressure_ = OneAtm;

    //list of reactors classes to choose from:
    //https://cantera.org/dev/cxx/db/d78/ReactorFactory_8cpp_source.html
    std::string reactor_type = "IdealGasReactor";

    //init times
    const double t0 = 0.0;
    const double tfinal = 1e-3;
    const double dt = 1e-6;
    const double dT0 = 1e-8; //1e-6 and lower

    //regular solution for jacobian
    auto sol = newSolution(mechanism_, name_, "none");
    auto gas = sol->thermo();
    gas->setState_TPX(temp_, pressure_, comp_mix);
    auto base_r = newReactor(reactor_type, sol);
    auto reactor = std::static_pointer_cast<Reactor>(base_r);
    ReactorNet net;
    net.addReactor(*reactor);

    //forward perturbation for central difference
    auto sol_f = newSolution(mechanism_, name_, "none");
    auto gas_f = sol_f->thermo();
    gas_f->setState_TPX(temp_ + dT0, pressure_, comp_mix);
    auto base_rf = newReactor(reactor_type, sol_f);
    auto r_f = std::static_pointer_cast<Reactor>(base_rf);
    ReactorNet net_f;
    net_f.addReactor(*r_f);

    //backward perturbation for central difference
    auto sol_b = newSolution(mechanism_, name_, "none");
    auto gas_b = sol_b->thermo();
    gas_b->setState_TPX(temp_ - dT0, pressure_, comp_mix);
    auto base_rb = newReactor(reactor_type, sol_b);
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
    std::cout << "no. of species = " << n_species << std::endl;
    fs::create_directories("bin");
    std::vector<std::ofstream> files(n_species);
    for (size_t k = 0; k < n_species; k++)
    {
        auto name = gas->speciesName(k);
        files[k].open("bin/" + name + ".csv");
        files[k]
            << "time,temp,"
            << "Y_" << name << ','
            << "eigen jacobian,"
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
        net.advance(t); //same as `integrator->integrate(tnow);`
        net_f.advance(t); //cent_diff
        net_b.advance(t); //cent_diff

        //get state solutions
        reactor->getState(y.data()); //same as `integrator->solution();`
        r_f->getState(yf.data());
        r_b->getState(yb.data());

        //calc eigen jacobian
        auto jac_start = std::chrono::high_resolution_clock::now(); //timing
        Eigen::SparseMatrix<double> eig_jac_sparse = reactor->finiteDifferenceJacobian();
        Eigen::MatrixXd eig_jac = Eigen::MatrixXd(eig_jac_sparse); //convert back to dense
        auto jac_end = std::chrono::high_resolution_clock::now(); //timing
        std::chrono::duration<double> jac_diff = jac_end - jac_start; //timting
        jac_total += jac_diff; //timing

        //SPECIES LOOP
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
                << t << std::scientific << ',' //time
                << y[0] << std::scientific << ',' //temp
                << Yk << std::scientific << ',' //species
                << jac_fd << std::scientific << ',' //jacobian eigen
                << central_diff << std::scientific << '\n'; //jacobian cent diff
        }
    }

    //print timing
    std::cout << "eigen jacobian time = " << jac_total.count() << " seconds" << std::endl;
    std::cout << "central difference time = " << cd_total.count() << " seconds" << std::endl;

    return 0;
}
