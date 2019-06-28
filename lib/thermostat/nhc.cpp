// Nose-Hoover Chain
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <chrono>

#include "../../include/simu_para.h"
#include "../../include/thermostat/nhc.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    void nhc_procedure::calc_cons_energy(void)
    {
        kine_energy = 0.0;
        pote_energy = 0.0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    kine_energy += pow(sys.molecules[mi].atoms[ai].p[di], 2)
                        / (2 * sys.molecules[mi].atoms[ai].m);
                    pote_energy += 0.5 * sys.molecules[mi].atoms[ai].m
                        * pow(omega, 2) * pow(sys.molecules[mi].atoms[ai].q[di], 2);
                }
            }
        }
        ther_energy = pow(tmvs[0].theta, 2) / (2 * tmvs[0].mu) + d * N * k * T * tmvs[0].eta;
        for (int j = 1; j < M - 1; j++) {
            ther_energy += pow(tmvs[j].theta, 2) / (2 * tmvs[j].mu) + k * T * tmvs[j].eta;
        }
        cons_energy = kine_energy + pote_energy + ther_energy;
    }

    // calculate physical forces (tmp, simple harmonic)
    void nhc_procedure::calc_physic_force(void)
    {
        for (auto mi = 0; mi < sys.molecules.size(); mi++)
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                for (auto di = 0; di < d; di++)
                    sys.molecules[mi].atoms[ai].F[di]
                    = -1 * sys.molecules[mi].atoms[ai].m
                    * pow(omega, 2) * sys.molecules[mi].atoms[ai].q[di];
    }

    void nhc_procedure::calc_thermo_force(int j)
    {
        if (j == 0) {
            kine_energy = 0.0;
            for (auto mi = 0; mi < sys.molecules.size(); mi++)
                for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                    for (auto di = 0; di < d; di++)
                        kine_energy += pow(sys.molecules[mi].atoms[ai].p[di], 2)
                        / (2 * sys.molecules[mi].atoms[ai].m);

            tmvs[0].Gamma = 2 * kine_energy - d * N * k * T;
        }
        else
            tmvs[j].Gamma = pow(tmvs[j - 1].theta, 2) / tmvs[j - 1].mu - k * T;
    }

    void nhc_procedure::physic_propagate(void)
    {
        calc_physic_force();
        for (auto mi = 0; mi < sys.molecules.size(); mi++)
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                for (auto di = 0; di < d; di++)
                    sys.molecules[mi].atoms[ai].p[di]
                    += bsp.time_step_size * sys.molecules[mi].atoms[ai].F[di] / 2;

        for (auto mi = 0; mi < sys.molecules.size(); mi++)
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                for (auto di = 0; di < d; di++)
                    sys.molecules[mi].atoms[ai].q[di]
                    += bsp.time_step_size * sys.molecules[mi].atoms[ai].p[di]
                    / sys.molecules[mi].atoms[ai].m;

        calc_physic_force();
        for (auto mi = 0; mi < sys.molecules.size(); mi++)
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                for (auto di = 0; di < d; di++)
                    sys.molecules[mi].atoms[ai].p[di]
                    += bsp.time_step_size * sys.molecules[mi].atoms[ai].F[di] / 2;
    }

    void nhc_procedure::thermo_propagate(void)
    {
        for (auto wi = 0; wi < tfs.n_sy; wi++) {
            double tmp_delta = tfs.weight[wi] * bsp.time_step_size / tfs.n_ff;
            for (auto ni = 0; ni < tfs.n_ff; ni++) {

                calc_thermo_force(M - 1);
                tmvs[M - 1].theta += tmp_delta * tmvs[M - 1].Gamma / 4;

                for (int j = M - 2; j >= 0; j--) {
                    tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                        / (8 * tmvs[j + 1].mu));
                    calc_thermo_force(j);
                    tmvs[j].theta += tmp_delta * tmvs[j].Gamma / 4;
                    tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                        / (8 * tmvs[j + 1].mu));
                }
                for (int j = 0; j < M; j++)
                    tmvs[j].eta += tmp_delta * tmvs[j].theta / (2 * tmvs[j].mu);

                double momen_scale = exp(-1 * tmp_delta * tmvs[0].theta / (2 * tmvs[0].mu));
                for (auto mi = 0; mi < sys.molecules.size(); mi++)
                    for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                        for (auto di = 0; di < d; di++)
                            sys.molecules[mi].atoms[ai].p[di] *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                        / (8 * tmvs[j + 1].mu));
                    calc_thermo_force(j);
                    tmvs[j].theta += tmp_delta * tmvs[j].Gamma / 4;
                    tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                        / (8 * tmvs[j + 1].mu));
                }
                calc_thermo_force(M - 1);
                tmvs[M - 1].theta += tmp_delta * tmvs[M - 1].Gamma / 4;
            }
        }
    }

    void nhc_procedure::print_nhc_procedure_title(std::ofstream& out) {
        std::cout << "\n\nNHC Procedure:\n   Time";
        out << "\n\nNHC Procedure:\nTime";

        for (int mi = 0; mi < sys.molecules.size(); mi++) {
            for (int ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (int di = 0; di < d; di++) {
                    std::cout << "    m[" << mi + 1 << "].a[" << ai + 1 << "].q[" << di + 1 << "]";
                    out << "\tm[" << mi + 1 << "].a[" << ai + 1 << "].q[" << di + 1 << "]";
                }
                for (int di = 0; di < d; di++) {
                    std::cout << "\tm[" << mi + 1 << "].a[" << ai + 1 << "].p[" << di + 1 << "]";
                    out << "\tm[" << mi + 1 << "].a[" << ai + 1 << "].p[" << di + 1 << "]";
                }
            }
        }
        std::cout << "\tcons_energy";
        out << "\tcons_energy";
    }

    void nhc_procedure::print_nhc_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        for (int mi = 0; mi < sys.molecules.size(); mi++) {
            for (int ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (int di = 0; di < d; di++) {
                    std::cout << std::setw(15) << sys.molecules[mi].atoms[ai].q[di];
                    out << "\t" << sys.molecules[mi].atoms[ai].q[di];
                }
                for (int di = 0; di < d; di++) {
                    std::cout << std::setw(15) << sys.molecules[mi].atoms[ai].p[di];
                    out << "\t" << sys.molecules[mi].atoms[ai].p[di];
                }
            }
        }
        std::cout << "\t" << cons_energy;
        out << "\t" << cons_energy;
    }
    
    void nhc_procedure::implement(std::ofstream& out, bool ifprint) {
        double t = 0;
        int ctr = 0;
        
        if (ifprint) {
            print_nhc_procedure_title(out);
            calc_cons_energy();
            print_nhc_procedure_data(out, t);
        }

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            // NHC numerical evolution
            thermo_propagate();
            physic_propagate();
            thermo_propagate();
            
            if (ifprint) {
                if (ctr == bsp.data_coll_peri) {
                    calc_cons_energy();
                    print_nhc_procedure_data(out, t);
                    ctr = 0;
                }
                ctr++;
            }
            t += bsp.time_step_size;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of NHC procedure (s): " << time_elapsed.count() << std::endl;

    }

}   // uhc
}   // thermostat
}   // uovie