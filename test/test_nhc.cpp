#include <string>
#include <vector>
#include <cmath>

#include <cassert>

#include "uov_proc.h"
#include "simu_para.h"
#include "phy_const.h"
#include "thermostat/nhc.h"

using namespace uovie;

uovie::Global::process simulation;

int main(int argc, char* argv[]) {

    /*** ================================================== ***/
    /*** Simulation Preparation                             ***/
    /*** ================================================== ***/

    simulation.open(argv[1]);
    simulation.read();

    /*** ================================================== ***/
    /*** Thermostat variables Generation                    ***/
    /*** ================================================== ***/

    std::vector<thermostat::nhc::thermo_vari> tmvs;
    thermostat::nhc::thermo_vari_generator(simulation.sys, tmvs, 4, 1);
    
    /*** ================================================== ***/
    /*** Thermostat Factorization Scheme Declaration        ***/
    /*** ================================================== ***/

    thermostat::nhc::thermo_factor_scheme tfs(7, 1);

    /*** ================================================== ***/
    /*** NHC Simulation Execution                           ***/
    /*** ================================================== ***/

    thermostat::nhc::nhc_procedure nhc_proce(simulation.bsp,
        simulation.sys, tmvs, tfs);
    nhc_proce.implement(simulation.out);

    /*** ================================================== ***/
    /*** Simulation termination                             ***/
    /*** ================================================== ***/
    
    simulation.print();
    simulation.close();

    return 0;
}