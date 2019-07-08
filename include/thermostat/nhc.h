/* Nose-Hoover Chain */
#ifndef NOSE_HOOVER_CHAIN_H_
#define NOSE_HOOVER_CHAIN_H_

// standard C++ headers
#include <fstream>
#include <vector>
#include <cmath>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "simu_para.h"
#include "phy_const.h"
#include "mol_geom.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    // shorthand for physical constants
    constexpr double h_bar = uovie::phy_const::red_Planck_const;
    constexpr double k = uovie::phy_const::Boltzmann_const;

    /*** ================================================== ***/
    /*** Thermostat Factorization Scheme                    ***/
    /*** ================================================== ***/

    // Suzuki¨CYoshida scheme
    class thermo_factor_scheme {
    public:
        thermo_factor_scheme() = default;
        thermo_factor_scheme(const int _n_sy, const int _n_ff): n_sy(_n_sy), n_ff(_n_ff) {
            if (_n_sy == 3)
                weight = { 1.351207191959658, -1.702414383919316, 1.351207191959658 };
            else if (_n_sy == 7)
                weight = { 0.784513610477560, 0.235573213359357, -1.17767998417887,
                    1.315186320683906, -1.17767998417887, 0.235573213359357, 0.784513610477560 };
            else
                throw "unsupported NHC factorization scheme";
        }

        int nsy() const { return n_sy; }
        int nff() const { return n_ff; }
        double w(const int& i) const { return weight[i]; }

    private:
        int n_sy;                       // the number of Suzuki-Yoshida weights
        int n_ff;                       // the number of terms in the further factorization
        std::vector<double> weight;     // Suzuki-Yoshida weights
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** NHC Procedure (Global)                             ***/
    /*** ================================================== ***/

    class nhc_procedure_global {
    public:
        nhc_procedure_global() = default;
        nhc_procedure_global(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const thermo_factor_scheme& _tfs, const int _nchain) :
            bsp(_bsp), sys(_sys), tfs(_tfs), nchain(_nchain) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const thermo_factor_scheme& tfs;
        const int nchain;   // extented dimension

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const int& M = nchain;

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force

        Eigen::ArrayXd mu;         // extented mass
        Eigen::ArrayXd eta;        // extented position
        Eigen::ArrayXd theta;      // extented momentum
        Eigen::ArrayXd Gamma;      // thermostat force

        double kine_energy = 0;
        double pote_energy = 0;
        double ther_energy = 0;
        double cons_energy = 0;

        void initialize();
        void calc_physic_force();
        void calc_thermo_force(const int& j);
        void physic_propagate();
        void thermo_propagate();
        void calc_cons_quant();

        void print_nhc_procedure_title(std::ofstream& out);
        void print_nhc_procedure_data(std::ofstream& out, double& t);
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** NHC Procedure (Local)                              ***/
    /*** ================================================== ***/

    class nhc_procedure_local {
    public:
        nhc_procedure_local() = default;
        nhc_procedure_local(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const thermo_factor_scheme& _tfs, const int _nchain) :
            bsp(_bsp), sys(_sys), tfs(_tfs), nchain(_nchain) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const thermo_factor_scheme& tfs;
        const int nchain;   // extented dimension

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const int& M = nchain;

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force

        Eigen::ArrayXXd mu;         // extented mass
        Eigen::ArrayXXd eta;        // extented position
        Eigen::ArrayXXd theta;      // extented momentum
        Eigen::ArrayXXd Gamma;      // thermostat force

        double kine_energy = 0;
        double pote_energy = 0;
        double ther_energy = 0;
        double cons_energy = 0;

        void initialize();
        void calc_physic_force();
        void calc_thermo_force(const int& j);
        void physic_propagate();
        void thermo_propagate();
        void calc_cons_quant();

        void print_nhc_procedure_title(std::ofstream& out);
        void print_nhc_procedure_data(std::ofstream& out, double& t);
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** NHC Procedure for pimd                             ***/
    /*** ================================================== ***/

    class nhc_procedure_for_pimd {
    public:
        nhc_procedure_for_pimd() = default;
        nhc_procedure_for_pimd(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const thermo_factor_scheme& _tfs, const int _nchain, const int _nbead):
            bsp(_bsp), sys(_sys), tfs(_tfs), nchain(_nchain), nbead(_nbead) { }
        
        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const thermo_factor_scheme& tfs;
        const int nchain;   // extented dimension
        const int nbead;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const int& M = nchain;
        const double fic_omega = sqrt(nbead) / (h_bar * beta);

        Eigen::ArrayXXd m;          // real mass
        Eigen::ArrayXXd q;          // real position
        Eigen::ArrayXXd m_tilde;    // transformed mass tilde
        Eigen::ArrayXXd m_bar;      // transformed mass bar
        Eigen::ArrayXXd r;          // transformed position
        Eigen::ArrayXXd s;          // fictitious momentum
        Eigen::ArrayXXd F;          // fictitious force
        
        Eigen::ArrayXXd mu;         // extented mass
        Eigen::ArrayXXd eta;        // extented position
        Eigen::ArrayXXd theta;      // extented momentum
        Eigen::ArrayXXd Gamma;      // thermostat force

        Eigen::ArrayXXd kine_energy;
        Eigen::ArrayXXd pote_energy;
        Eigen::ArrayXXd ther_energy;
        Eigen::ArrayXXd cons_energy;

        double prim_kine_estor = 0;
        double prim_pote_estor = 0;
        double prim_pres_estor = 0;

        void initialize();
        void stag_trans();
        void inve_stag_trans();
        void calc_physic_force();
        void calc_thermo_force(const int& j);
        void physic_propagate();
        void thermo_propagate();
        void calc_cons_quant();
        void calc_prim_estor();

        void print_nhc_procedure_title(std::ofstream& out);
        void print_nhc_procedure_data(std::ofstream& out, double& t);
    };

} // !nhc
} // !thermostat
} // !uovie
#endif // !NOSE_HOOVER_CHAIN_H_