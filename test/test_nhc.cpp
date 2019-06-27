#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "../include/iostream.h"
#include "../include/simu_para.h"
#include "../include/phy_const.h"
#include "../include/thermostat/nose_hoover_chain.h"

using namespace uovie;

uovie::Global::process cano_ens_nhc_simu;

int main(int argc, char* argv[]) {

    /*** ================================================== ***/
    /*** Prepare a Canonical Ensemble NHC Simulation        ***/
    /*** ================================================== ***/

    std::string filename = argv[1];
    cano_ens_nhc_simu.open(filename);
    cano_ens_nhc_simu.read();

    double T = 10 / Phy_Const::Boltzmann_const; // temperature
    double omega = 1;

    int d = cano_ens_nhc_simu.sys.dimension;

    int N = 0;                              // system dimension
    for (auto i = 0; i < cano_ens_nhc_simu.sys.molecules.size(); i++)
        N += static_cast<int>(cano_ens_nhc_simu.sys.molecules[i].atoms.size());

    /*** ================================================== ***/
    /*** Declare Basic Simulation Parameters                ***/
    /*** ================================================== ***/

    Global::basic_simu_para bsp(5000.0, 0.0001, 1000);

    /*** ================================================== ***/
    /*** Declare Thermostat Variables                       ***/
    /*** ================================================== ***/

    const int M = 4;                        // extented dimension
    std::vector<thermostat::nhc::thermo_vari> tmvs;

    tmvs.push_back({ static_cast<double>(d * N), T, 0, 1 });
    tmvs.push_back({ 1, T, 0, -1 });
    tmvs.push_back({ 1, T, 0, 1 });
    tmvs.push_back({ 1, T, 0, -1 });
    //for (int j = 1; j < M; j++)
    //    tmvs.push_back({ 1, T, 0, 1 });

    /*** ================================================== ***/
    /*** Declare a Suzuki¨CYoshida scheme and initialize it ***/
    /*** ================================================== ***/

    thermostat::nhc::thermo_factor_scheme six_ord_sy_scheme(7, 1);
    six_ord_sy_scheme.init();

    /*** ================================================== ***/
    /*** Execute an NHC Simulation                          ***/
    /*** ================================================== ***/

    thermostat::nhc::nhc_procedure execute_nhc_proce(cano_ens_nhc_simu.sys,
        T, tmvs, six_ord_sy_scheme, bsp);
    execute_nhc_proce.implement(cano_ens_nhc_simu.out);

    cano_ens_nhc_simu.print();
    cano_ens_nhc_simu.close();

    return 0;
}