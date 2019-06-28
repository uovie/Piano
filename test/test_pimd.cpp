#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "../include/iostream.h"
#include "../include/simu_para.h"
#include "../include/phy_const.h"
#include "../include/pimd.h"

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
    /*** PIMD Simulation Execution                          ***/
    /*** ================================================== ***/

    pimd::nhc_procedure pimd_proce(simulation.sys,
        T, tmvs, six_ord_sy_scheme, bsp);
    pimd_proce.implement(simulation.out);

    simulation.print();
    simulation.close();

    return 0;

}