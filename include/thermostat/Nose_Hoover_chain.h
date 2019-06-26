/* Nose-Hoover Chain */
#ifndef NOSE_HOOVER_CHAIN_H_
#define NOSE_HOOVER_CHAIN_H_

// standard C++ headers
#include <fstream>
#include <vector>
#include <cmath>

// uovie headers
#include "../phy_const.h"
#include "../mol_geom.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    /*** ==================== ***/
    /*** Thermostat Variables ***/
    /*** ==================== ***/

    class thermo_vari {
    public:
        thermo_vari() = default;
        thermo_vari(double dof, double T, double e, double t):
            mu(dof* Phy_Const::Boltzmann_const* T), eta(e), theta(t) { }

        double mu;                      // extented mass
        double eta;                     // extented position
        double theta;                   // extented momentum
        double Gamma;                   // thermostat force
    };
    
    /*** =============================== ***/
    /*** Thermostat Factorization Scheme ***/
    /*** =============================== ***/

    // Suzuki¨CYoshida scheme
    class thermo_factor_scheme {
    public:
        thermo_factor_scheme() = default;
        thermo_factor_scheme(const int nsy, const int nff): n_sy(nsy), n_ff(nff) { }

        int n_sy;                       // the number of Suzuki-Yoshida weights
        int n_ff;                       // the number of terms in the further factorization
        std::vector<double> weight;     // Suzuki-Yoshida weights

        void init() {
            if (n_sy == 3)
                weight = { 1.351207191959658, -1.702414383919316, 1.351207191959658 };
            else if (n_sy == 7)
                weight = { 0.784513610477560, 0.235573213359357, -1.17767998417887,
                    1.315186320683906, -1.17767998417887, 0.235573213359357, 0.784513610477560 };
            else
                throw "unsupported Suzuki-Yoshida scheme";
        }
    };

    /*** ==================== ***/
    /*** NHC Procedure        ***/
    /*** ==================== ***/

    class nhc_procedure {
    public:
        nhc_procedure() = default;
        nhc_procedure(Global::system& s, double& te, std::vector<thermo_vari>& tvs,
            thermo_factor_scheme& fs, Global::basic_simu_para& sp):
            sys(s), T(te), tmvs(tvs), tfs(fs), bsp(sp) { }

        Global::system sys;
        double T;                                       // temperature
        std::vector<thermo_vari> tmvs;
        thermo_factor_scheme tfs;
        Global::basic_simu_para bsp;
        
        void implement(std::ofstream& out);

    private:
        int N = 0;                                      // system dimension
        void init() {
            for (auto i = 0; i < sys.molecules.size(); i++)
                N += static_cast<int>(sys.molecules[i].atoms.size());
        }
        int M = static_cast<int>(tmvs.size());          // extented dimension

        // calculate physical forces (tmp, simple harmonic)
        void physic_force(void)
        {
            double omega = 1;

            for (auto mi = 0; mi < sys.molecules.size(); mi++)
                for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                    for (auto di = 0; di < 3; di++)
                        sys.molecules[mi].atoms[ai].F[di]
                            = -1 * sys.molecules[mi].atoms[ai].m
                            * pow(omega, 2) * sys.molecules[mi].atoms[ai].q[di];
        }

        // calculate thermostat forces
        void thermo_force(int j)
        {
            if (j == 0) {
                double kinetic_energy = 0.0;

                for (auto mi = 0; mi < sys.molecules.size(); mi++)
                    for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                        for (auto di = 0; di < 3; di++)
                            kinetic_energy += pow(sys.molecules[mi].atoms[ai].p[di], 2)
                            / (2 * sys.molecules[mi].atoms[ai].m);

                tmvs[0].Gamma = 2 * kinetic_energy - 3 * N * Phy_Const::Boltzmann_const * T;
            }
            else
                tmvs[j].Gamma = pow(tmvs[j - 1].theta, 2) / tmvs[j - 1].mu - Phy_Const::Boltzmann_const * T;
        }

        // physical propagation
        void physic_propagate(void)
        {
            physic_force();
            for (auto mi = 0; mi < sys.molecules.size(); mi++)
                for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                    for (auto di = 0; di < 3; di++)
                        sys.molecules[mi].atoms[ai].p[di]
                        += bsp.step_size * sys.molecules[mi].atoms[ai].F[di] / 2;

            for (auto mi = 0; mi < sys.molecules.size(); mi++)
                for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                    for (auto di = 0; di < 3; di++)
                        sys.molecules[mi].atoms[ai].q[di]
                        += bsp.step_size * sys.molecules[mi].atoms[ai].p[di]
                        / sys.molecules[mi].atoms[ai].m;

            physic_force();
            for (auto mi = 0; mi < sys.molecules.size(); mi++)
                for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                    for (auto di = 0; di < 3; di++)
                        sys.molecules[mi].atoms[ai].p[di]
                        += bsp.step_size * sys.molecules[mi].atoms[ai].F[di] / 2;
        }

        // thermostat propagation
        void thermo_propagate(void)
        {
            for (auto wi = 0; wi < tfs.n_sy; wi++) {
                double tmp_delta = tfs.weight[wi] * bsp.step_size / tfs.n_ff;
                for (auto ni = 0; ni < tfs.n_ff; ni++) {

                    thermo_force(M - 1);
                    tmvs[M - 1].theta += tmp_delta * tmvs[M - 1].Gamma / 4;

                    for (int j = M - 2; j >= 0; j--) {
                        tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                            / (8 * tmvs[j + 1].mu));
                        thermo_force(j);
                        tmvs[j].theta += tmp_delta * tmvs[j].Gamma / 4;
                        tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                            / (8 * tmvs[j + 1].mu));
                    }
                    for (int j = 0; j < M; j++)
                        tmvs[j].eta += tmp_delta * tmvs[j].theta / (2 * tmvs[j].mu);

                    double momen_scale = exp(-1 * tmp_delta * tmvs[0].theta / (2 * tmvs[0].mu));
                    for (auto mi = 0; mi < sys.molecules.size(); mi++)
                        for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++)
                            for (auto di = 0; di < 3; di++)
                                sys.molecules[mi].atoms[ai].p[di] *= momen_scale;
                    
                    for (int j = 0; j < M - 1; j++) {
                        tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                            / (8 * tmvs[j + 1].mu));
                        thermo_force(j);
                        tmvs[j].theta += tmp_delta * tmvs[j].Gamma / 4;
                        tmvs[j].theta *= exp(-1 * tmp_delta * tmvs[j + 1].theta
                            / (8 * tmvs[j + 1].mu));
                    }
                    thermo_force(M-1);
                    tmvs[M - 1].theta += tmp_delta * tmvs[M - 1].Gamma / 4;
                }
            }
        }

    };
        
}   // nhc
}   // thermostat
}   // uovie

#endif