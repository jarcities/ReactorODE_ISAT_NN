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
