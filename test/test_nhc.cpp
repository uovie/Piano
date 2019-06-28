#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "../include/process.h"
#include "../include/simu_para.h"
#include "../include/phy_const.h"
#include "../include/thermostat/nhc.h"

using namespace uovie;

uovie::Global::process simulation;

int main(int argc, char* argv[]) {

    /*** ================================================== ***/
    /*** Simulation Preparation                             ***/
    /*** ================================================== ***/

    simulation.open(argv[1]);
    simulation.read();

    int& d = simulation.sys.dimension;
    int& N = simulation.sys.num_part;
    double& T = simulation.sys.temperature;

    /*** ================================================== ***/
    /*** Thermostat variables Declaration                   ***/
    /*** ================================================== ***/

    std::vector<thermostat::nhc::thermo_vari> tmvs;
    tmvs.push_back({ static_cast<double>(d * N), T, 0, 1 });
    //tmvs.push_back({ 1, T, 0, -1 });
    //tmvs.push_back({ 1, T, 0, 1 });
    //tmvs.push_back({ 1, T, 0, -1 });

    /*** ================================================== ***/
    /*** Thermostat Factorization Scheme Declaration        ***/
    /*** ================================================== ***/

    thermostat::nhc::thermo_factor_scheme tfs(7, 1);
    tfs.init();

    /*** ================================================== ***/
    /*** NHC Simulation Execution                           ***/
    /*** ================================================== ***/

    thermostat::nhc::nhc_procedure nhc_proce(simulation.bsp,
        simulation.sys, tmvs, tfs);
    nhc_proce.implement(simulation.out, 1);

    /*** ================================================== ***/
    /*** Simulation termination                             ***/
    /*** ================================================== ***/
    
    simulation.print();
    simulation.close();

    return 0;
}