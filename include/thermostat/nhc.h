/* Nose-Hoover Chain */
#ifndef NOSE_HOOVER_CHAIN_H_
#define NOSE_HOOVER_CHAIN_H_

// standard C++ headers
#include <fstream>
#include <vector>
#include <cmath>

// uovie headers
#include "../simu_para.h"
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
            mu(dof* phy_const::Boltzmann_const* T), eta(e), theta(t) { }

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
        nhc_procedure(Global::basic_simu_para& b, Global::system& s,
            std::vector<thermo_vari>& tvs, thermo_factor_scheme& fs) :
            bsp(b), sys(s), tmvs(tvs), tfs(fs) { }

        Global::basic_simu_para bsp;
        Global::system sys;
        std::vector<thermo_vari> tmvs;
        thermo_factor_scheme tfs;

        void implement(std::ofstream& out, bool ifprint);

    private:
        int& d = sys.dimension;
        int& N = sys.num_part;
        double k = phy_const::Boltzmann_const;
        double& T = sys.temperature;
        int M = tmvs.size();          // extented dimension

        double kine_energy;
        double pote_energy;
        double ther_energy;
        double cons_energy;

        double omega = 1;

        void calc_cons_energy(void);
        void calc_physic_force(void);
        void calc_thermo_force(int j);
        void physic_propagate(void);
        void thermo_propagate(void);

        void print_nhc_procedure_title(std::ofstream& out);
        void print_nhc_procedure_data(std::ofstream& out, double& t);

    };
}   // nhc
}   // thermostat
}   // uovie
#endif