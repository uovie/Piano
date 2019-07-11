/* Path Integral Molecular Dynamics */
#ifndef PIMD_H_
#define PIMD_H_

// standard C++ headers
#include <fstream>
#include <random>
#include <chrono>

// uovie headers
#include "simu_para.h"
#include "mol_geom.h"
#include "phy_const.h"
#include "thermostat/nhc.h"

namespace uovie {
namespace pimd {

    // shorthand for physical constants
    constexpr double h_bar = uovie::phy_const::red_Planck_const;
    constexpr double k = uovie::phy_const::Boltzmann_const;
    
    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** PIMD procedure (Base)                              ***/
    /*** ================================================== ***/

    class pimd_base {
    public:
        pimd_base() = default;
        pimd_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead) :
            fn_no_ex(_fn_no_ex), bsp(_bsp), sys(_sys), nbead(_nbead) { }

        void implement();

    protected:
        const std::string& fn_no_ex;
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const int& nbead;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double fic_omega = sqrt(nbead) / (h_bar * beta);

        Eigen::ArrayXXd m;          // real mass
        Eigen::ArrayXXd q;          // real position
        Eigen::ArrayXXd m_tilde;    // transformed mass tilde
        Eigen::ArrayXXd m_bar;      // transformed mass bar
        Eigen::ArrayXXd r;          // transformed position
        Eigen::ArrayXXd s;          // fictitious momentum
        Eigen::ArrayXXd F;          // fictitious force

        Eigen::ArrayXXd kin_ene;
        Eigen::ArrayXXd pot_ene;

        double prim_kin_estor = 0;
        double prim_pot_estor = 0;
        double prim_pre_estor = 0;

        virtual void initialize() = 0;
        void stag_trans();
        void inve_stag_trans();
        void calc_physic_force();
        void calc_prim_estor();
        virtual void implement_one_step() = 0;

        std::vector<std::string> labels;

        virtual void print_pimd_proce_title(std::ofstream& chk, std::ofstream& out);
        virtual void print_pimd_proce_data(std::ofstream& chk, std::ofstream& out, double& t);
        virtual void print_conclusion_info(std::ofstream& chk, std::ofstream& out,
            const std::chrono::duration<double>& time_elap);
    };
    
    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** PIMD procedure via LD (Base)                     ***/
    /*** ================================================== ***/

    class pimd_via_ld_base : public pimd_base {
    public:
        pimd_via_ld_base() = default;
        pimd_via_ld_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const double& _gamma):
            pimd_base(_fn_no_ex, _bsp, _sys, _nbead), gamma(_gamma) { }

    protected:
        const double& gamma;    // friction coefficient

        Eigen::ArrayXXd nrand;
        std::vector<std::mt19937> mtes;
        std::vector<std::normal_distribution<double>> nds;

        void initialize() override;
        void upd_rand_num();
    };

    /*** ================================================== ***/
    /*** PIMD procedure via LD (Side)                       ***/
    /*** ================================================== ***/

    class pimd_via_ld_side : public pimd_via_ld_base {
    public:
        pimd_via_ld_side() = default;
        pimd_via_ld_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const double& _gamma) :
            pimd_via_ld_base(_fn_no_ex, _bsp, _sys, _nbead, _gamma)
        {
            labels.push_back("LD");
            labels.push_back("Side");
        }

    private:
        const double c1 = exp(-gamma * Dt / 2);
        const double c2 = sqrt(1 - pow(c1, 2));

        void implement_one_step() override;
    };

    /*** ================================================== ***/
    /*** PIMD procedure via LD (Middle)                     ***/
    /*** ================================================== ***/

    class pimd_via_ld_middle : public pimd_via_ld_base {
    public:
        pimd_via_ld_middle() = default;
        pimd_via_ld_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const double& _gamma) :
            pimd_via_ld_base(_fn_no_ex, _bsp, _sys, _nbead, _gamma)
        {
            labels.push_back("LD");
            labels.push_back("Middle");
        }

    private:
        const double c1 = exp(-gamma * Dt);
        const double c2 = sqrt(1 - pow(c1, 2));

        void implement_one_step() override;
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** PIMD procedure via AT (Base)                       ***/
    /*** ================================================== ***/

    class pimd_via_at_base : public pimd_base {
    public:
        pimd_via_at_base() = default;
        pimd_via_at_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const double& _nu) :
            pimd_base(_fn_no_ex, _bsp, _sys, _nbead), nu(_nu) { }

    protected:
        const double& nu;

        Eigen::ArrayXXd nrand;
        Eigen::ArrayXXd urand;

        std::vector<std::mt19937> nmtes;
        std::vector<std::normal_distribution<double>> nds;
        std::vector<std::mt19937> umtes;
        std::vector<std::uniform_real_distribution<double>> urds;

        void initialize() override;
        void upd_uni_rand_num();
    };

    /*** ================================================== ***/
    /*** PIMD procedure via AT (Side)                       ***/
    /*** ================================================== ***/

    class pimd_via_at_side : public pimd_via_at_base {
    public:
        pimd_via_at_side() = default;
        pimd_via_at_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const double& _nu) :
            pimd_via_at_base(_fn_no_ex, _bsp, _sys, _nbead, _nu)
        {
            labels.push_back("AT");
            labels.push_back("Side");
        }

    private:
        const double cri = 1 - exp(-nu * Dt / 2);
        void implement_one_step();
    };

    /*** ================================================== ***/
    /*** PIMD procedure via AT (Middle)                     ***/
    /*** ================================================== ***/

    class pimd_via_at_middle : public pimd_via_at_base {
    public:
        pimd_via_at_middle() = default;
        pimd_via_at_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const double& _nu) :
            pimd_via_at_base(_fn_no_ex, _bsp, _sys, _nbead, _nu)
        {
            labels.push_back("AT");
            labels.push_back("Middle");
        }

    private:
        const double cri = 1 - exp(-nu * Dt);
        void implement_one_step();
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** PIMD via NHC (Base)                                ***/
    /*** ================================================== ***/

    class pimd_via_nhc_base : public pimd_base {
    public:
        pimd_via_nhc_base() = default;
        pimd_via_nhc_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int& _nbead, const thermostat::nhc::thermo_factor_scheme& _tfs,
            const int _nchain) : pimd_base(_fn_no_ex, _bsp, _sys, _nbead), tfs(_tfs), nchain(_nchain) { }

    protected:
        const thermostat::nhc::thermo_factor_scheme& tfs;
        const int nchain;   // extented dimension
        const int& M = nchain;

        Eigen::ArrayXXd mu;         // extented mass
        Eigen::ArrayXXd eta;        // extented position
        Eigen::ArrayXXd theta;      // extented momentum
        Eigen::ArrayXXd Gamma;      // thermostat force

        Eigen::ArrayXXd kin_ene_arr;
        Eigen::ArrayXXd pot_ene_arr;
        Eigen::ArrayXXd the_ene_arr;
        Eigen::ArrayXXd con_ene_arr;

        void initialize() override;
        void calc_thermo_force(const int& j);
        void calc_cons_quant();

        void print_pimd_proce_title(std::ofstream& chk, std::ofstream& out) override;
        void print_pimd_proce_data(std::ofstream& chk, std::ofstream& out, double& t) override;
    };

    /*** ================================================== ***/
    /*** PIMD via NHC (Side)                                ***/
    /*** ================================================== ***/

    class pimd_via_nhc_side : public pimd_via_nhc_base {
    public:
        pimd_via_nhc_side() = default;
        pimd_via_nhc_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int _nbead, const thermostat::nhc::thermo_factor_scheme& _tfs,
            const int _nchain) : pimd_via_nhc_base(_fn_no_ex, _bsp, _sys, _nbead, _tfs, _nchain)
        {
            labels.push_back("NHC");
            labels.push_back("Side");
        }

    private:
        void implement_one_step();
    };

    /*** ================================================== ***/
    /*** PIMD via NHC (Middle)                                ***/
    /*** ================================================== ***/

    class pimd_via_nhc_middle : public pimd_via_nhc_base {
    public:
        pimd_via_nhc_middle() = default;
        pimd_via_nhc_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int _nbead, const thermostat::nhc::thermo_factor_scheme& _tfs,
            const int _nchain) : pimd_via_nhc_base(_fn_no_ex, _bsp, _sys, _nbead, _tfs, _nchain)
        {
            labels.push_back("NHC");
            labels.push_back("Middle");
        }

    private:
        void implement_one_step();
    };

} // !pimd
} // !uovie
#endif // !PIMD_H_