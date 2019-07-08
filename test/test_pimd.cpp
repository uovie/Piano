#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "uov_proc.h"
#include "simu_para.h"
#include "phy_const.h"
#include "pimd.h"

using namespace uovie;

uovie::Global::process simulation;

int main(int argc, char* argv[]) {

    /*** ================================================== ***/
    /*** Simulation Preparation                             ***/
    /*** ================================================== ***/

    simulation.open(argv[1]);
    simulation.read();

    /*** ================================================== ***/
    /*** PIMD Simulation Execution                          ***/
    /*** ================================================== ***/

    pimd::pimd_procedure pimd_proce(simulation.bsp, simulation.sys, 6);
    pimd_proce.implement(simulation.out, simulation.fn_no_ex);

    /*** ================================================== ***/
    /*** Simulation termination                             ***/
    /*** ================================================== ***/

    simulation.print();
    simulation.close();

    return 0;

}