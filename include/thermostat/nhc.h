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
        thermo_vari(double m, double e, double t):
            mu(m), eta(e), theta(t) { }

        double mu;                      // extented mass
        double eta;                     // extented position
        double theta;                   // extented momentum
        double Gamma;                   // thermostat force
    };

    void thermo_vari_generator(const Global::system& sys,
        std::vector<thermo_vari>& tmvs, int M, double tau);
    
    /*** =============================== ***/
    /*** Thermostat Factorization Scheme ***/
    /*** =============================== ***/

    // Suzuki¨CYoshida scheme
    class thermo_factor_scheme {
    public:
        thermo_factor_scheme() = default;
        thermo_factor_scheme(const int nsy, const int nff): n_sy(nsy), n_ff(nff) {
            if (nsy == 3)
                weight = { 1.351207191959658, -1.702414383919316, 1.351207191959658 };
            else if (nsy == 7)
                weight = { 0.784513610477560, 0.235573213359357, -1.17767998417887,
                    1.315186320683906, -1.17767998417887, 0.235573213359357, 0.784513610477560 };
            else
                throw "unsupported Suzuki-Yoshida scheme";
        }

        int nsy() const { return n_sy; }
        int nff() const { return n_ff; }
        double w(const int& i) const { return weight[i]; }

    private:
        int n_sy;                       // the number of Suzuki-Yoshida weights
        int n_ff;                       // the number of terms in the further factorization
        std::vector<double> weight;     // Suzuki-Yoshida weights
    };

    /*** ================================================== ***/
    /*** NHC Procedure Base                                 ***/
    /*** ================================================== ***/

    class nhc_procedure_base {
    public:
        nhc_procedure_base() = default;
        nhc_procedure_base(const Global::basic_simu_para& b, Global::system& s,
            std::vector<thermo_vari>& tvs, const thermo_factor_scheme& fs) :
            bsp(b), sys(s), tmvs(tvs), tfs(fs) { }
        //~nhc_procedure_base();
        
        const Global::basic_simu_para& bsp;
        Global::system& sys;
        std::vector<thermo_vari>& tmvs;
        const thermo_factor_scheme& tfs;
        
        void implement();
        void implement(std::ofstream& out);

    private:
        double omega = 1;

    protected:
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const double k = phy_const::Boltzmann_const;
        const double& T = sys.temperature;
        const int M = tmvs.size();          // extented dimension

        double kine_energy;
        double pote_energy;
        double ther_energy;
        double cons_energy;

        virtual void calc_physic_force();
        void calc_thermo_force(const int& j);
        void physic_propagate();
        void thermo_propagate();

        virtual void calc_cons_energy();

        void print_nhc_procedure_title(std::ofstream& out);
        void print_nhc_procedure_data(std::ofstream& out, double& t);
    };

    /*** ================================================== ***/
    /*** NHC Procedure for PIMD                             ***/
    /*** ================================================== ***/

    //class nhc_procedure_for_pimd : public nhc_procedure_base {
    //public:
    //    nhc_procedure_for_pimd() = default;
    //    nhc_procedure_for_pimd(Global::basic_simu_para& b, Global::system& s,
    //        std::vector<thermo_vari>& tvs, thermo_factor_scheme& fs, const int& i) :
    //        bsp(b), sys(s), tmvs(tvs), tfs(fs), bi(i) { }

    //private:
    //    const int bi;    // bead index
    //    double fic_omega;

    //    void calc_physic_force() override;
    //    void calc_cons_energy() override;
    //};

}   // nhc
}   // thermostat
}   // uovie
#endif