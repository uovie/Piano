#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>

#include "pimd.h"

namespace uovie {
namespace pimd {

    /*** ================================================== ***/
    /*** Staging Transformation                             ***/
    /*** ================================================== ***/

    std::vector<Global::system> prim_sys_coll::stag_trans()
    {
        std::vector<Global::system> stsyss;
        stsyss.push_back(syss[0]);
        for (auto bi = 1; bi < syss.size() - 1; bi++) {
            stsyss.push_back(syss[0]);
            for (auto mi = 0; mi < syss[bi].molecules.size(); mi++) {
                for (auto ai = 0; ai < syss[bi].molecules[mi].atoms.size(); ai++) {
                    stsyss[bi].molecules[mi].atoms[ai].m
                        = ((bi + 1) / bi) * syss[bi].molecules[mi].atoms[ai].m;
                    for (auto di = 0; di < syss[bi].dimension; di++) {
                        stsyss[bi].molecules[mi].atoms[ai].q[di]
                            = -(1 / (bi + 1)) * syss[0].molecules[mi].atoms[ai].q[di]
                            + syss[bi].molecules[mi].atoms[ai].q[di]
                            - (bi / (bi + 1)) * syss[bi + 1].molecules[mi].atoms[ai].q[di];
                    }
                }
            }
        }

        stsyss.push_back(syss[0]);
        int last_bi = syss.size() - 1;
        for (auto mi = 0; mi < syss[last_bi].molecules.size(); mi++) {
            for (auto ai = 0; ai < syss[last_bi].molecules[mi].atoms.size(); ai++) {
                stsyss[last_bi].molecules[mi].atoms[ai].m
                    = ((last_bi + 1) / last_bi) * syss[last_bi].molecules[mi].atoms[ai].m;
                for (auto di = 0; di < syss[last_bi].dimension; di++) {
                    stsyss[last_bi].molecules[mi].atoms[ai].q[di]
                        = syss[last_bi].molecules[mi].atoms[ai].q[di]
                        - syss[0].molecules[mi].atoms[ai].q[di];
                }
            }
        }

        return stsyss;
    }

    /*** ================================================== ***/
    /*** Inverse Staging Transformation                     ***/
    /*** ================================================== ***/

    std::vector<Global::system> stra_sys_coll::inve_stag_trans()
    {
        std::vector<Global::system> prsyss;
        prsyss.push_back(syss[0]);
        for (auto bi = 1; bi < syss.size(); bi++) {
            prsyss.push_back(syss[0]);
            for (auto mi = 0; mi < syss[bi].molecules.size(); mi++) {
                for (auto ai = 0; ai < syss[bi].molecules[mi].atoms.size(); ai++) {
                    for (auto di = 0; di < syss[bi].dimension; di++) {
                        prsyss[bi].molecules[mi].atoms[ai].q[di]
                            = syss[0].molecules[mi].atoms[ai].q[di];
                        for (auto bj = bi; bj < syss.size(); bj++)
                            prsyss[bi].molecules[mi].atoms[ai].q[di]
                                += (bi / bj) * syss[bj].molecules[mi].atoms[ai].q[di];
                    }
                }
            }
        }
        return prsyss;
    }

    void pimd_procedure::implement()
    {
        /* NHC Preparation */
        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        std::vector<std::vector<thermostat::nhc::thermo_vari>> tmvsc(nbead);
        std::vector<thermostat::nhc::nhc_procedure_for_pimd> nhc_proces;
        for (int bi = 0; bi < nbead; bi++) {
            thermostat::nhc::thermo_vari_generator(stsc.syss[bi], tmvsc[bi], 4, 1);
            nhc_proces.push_back(thermostat::nhc::nhc_procedure_for_pimd(
                bsp, sys, tmvsc[bi], tfs, prsc.syss, stsc.syss, bi));
        }
        std::cout << "NHC is ready, and PIMD starts." << std::endl;
        /* NHC Preparation */

        for (double t = 0; t <= bsp.run_time; t += bsp.time_step_size) {
            for (int bi = 0; bi < nbead; bi++) {
                nhc_proces[bi].imple_one_step();
                prsc.syss = stsc.inve_stag_trans();
            }
            t += bsp.time_step_size;
        }
    }

    void pimd_procedure::implement(std::ofstream& out, std::string& fn_no_ex)
    {
        double t = 0;
        int ctr = 0;

        /* NHC Preparation */
        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        std::vector<thermostat::nhc::thermo_vari_coll> tmvcs(nbead);
        std::vector<thermostat::nhc::nhc_procedure_for_pimd> nhc_proces;
        for (int bi = 0; bi < nbead; bi++) {
            thermostat::nhc::thermo_vari_generator(stsc.syss[bi], tmvcs[bi].tmvc, 4, 1);
            nhc_proces.push_back(thermostat::nhc::nhc_procedure_for_pimd(
                bsp, sys, tmvcs[bi].tmvc, tfs, prsc.syss, stsc.syss, bi));
        }
        std::cout << "NHC is ready, and PIMD starts." << std::endl;
        /* NHC Preparation */

        std::ofstream out_cpd_data;
        out_cpd_data.open(fn_no_ex + ".cpd.dat");

        print_pimd_cpd_title(out_cpd_data);
        tot_sys_energy = 0;
        for (int bi = 0; bi < nbead; bi++) {
            nhc_proces[bi].calc_syco_energy();
            tot_sys_energy += nhc_proces[bi].sys_ene();
        }
        cano_prob_dens = exp(-tot_sys_energy / (k * T));
        print_pimd_cpd_data(out_cpd_data, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            for (int bi = 0; bi < nbead; bi++)
                nhc_proces[bi].imple_one_step();

            prsc.syss = stsc.inve_stag_trans();

            if (ctr == bsp.data_coll_peri) {
                tot_sys_energy = 0;
                for (int bi = 0; bi < nbead; bi++) {
                    nhc_proces[bi].calc_syco_energy();
                    tot_sys_energy += nhc_proces[bi].sys_ene();
                }
                cano_prob_dens = exp(-tot_sys_energy / (k * T));
                print_pimd_cpd_data(out_cpd_data, t);
                ctr = 0;
            }
            
            ctr++;
            t += bsp.time_step_size;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;
        std::cout << "\nElapsed time of PIMD procedure (s): " << time_elapsed.count() << std::endl;
    }

    void pimd_procedure::print_pimd_cpd_title(std::ofstream& out) {
        std::cout << "\n\nPIMD Procedure:\n   Time";
        out << "\n\nPIMD Procedure:\nTime";
        std::cout << "    cano_prob_dens";
        out << "\tcano_prob_dens";
    }

    void pimd_procedure::print_pimd_cpd_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::setw(15) << cano_prob_dens;
        out << "\t" << cano_prob_dens;
    }

} // !pimd
} // !uovie