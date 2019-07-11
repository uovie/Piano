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
#include "thermostat/thermostat.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    //--------------------------------------------------------//

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
    /*** NHC Procedure (Base)                               ***/
    /*** ================================================== ***/

    class nhc_base : public thermostat_base {
    public:
        nhc_base() = default;
        nhc_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            thermostat_base(_fn_no_ex, _bsp, _sys), tfs(_tfs), nchain(_nchain) { }

    protected:
        const thermo_factor_scheme& tfs;
        const int nchain;   // extented dimension
        const int& M = nchain;

        double the_ene = 0;
        double con_ene = 0;

        void print_ther_proce_title(std::ofstream& chk, std::ofstream& out) override;
        void print_ther_proce_data(std::ofstream& chk, std::ofstream& out, double& t) override;
        void print_conclusion_info(std::ofstream& chk, std::ofstream& out,
            const std::chrono::duration<double>& time_elap) override;
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** NHC Procedure (Global, Base)                       ***/
    /*** ================================================== ***/

    class nhc_global_base : public nhc_base {
    public:
        nhc_global_base() = default;
        nhc_global_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            nhc_base(_fn_no_ex, _bsp, _sys, _tfs, _nchain) { }

    protected:
        Eigen::ArrayXd mu;         // extented mass
        Eigen::ArrayXd eta;        // extented position
        Eigen::ArrayXd theta;      // extented momentum
        Eigen::ArrayXd Gamma;      // thermostat force

        void initialize() override;
        void calc_thermo_force(const int& j);
        void calc_cons_quant();
    };

    /*** ================================================== ***/
    /*** NHC Procedure (Global, Side)                       ***/
    /*** ================================================== ***/

    class nhc_global_side : public nhc_global_base {
    public:
        nhc_global_side() = default;
        nhc_global_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            nhc_global_base(_fn_no_ex, _bsp, _sys, _tfs, _nchain)
        {
            labels.push_back("NHC");
            labels.push_back("Global");
            labels.push_back("Side");
        }

    private:
        void implement_one_step() override;
    };

    /*** ================================================== ***/
    /*** NHC Procedure (Global, Middle)                     ***/
    /*** ================================================== ***/

    class nhc_global_middle : public nhc_global_base {
    public:
        nhc_global_middle() = default;
        nhc_global_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            nhc_global_base(_fn_no_ex, _bsp, _sys, _tfs, _nchain)
        {
            labels.push_back("NHC");
            labels.push_back("Global");
            labels.push_back("Middle");
        }

    private:
        void implement_one_step() override;
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** NHC Procedure (Local, Base)                        ***/
    /*** ================================================== ***/

    class nhc_local_base : public nhc_base {
    public:
        nhc_local_base() = default;
        nhc_local_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            nhc_base(_fn_no_ex, _bsp, _sys, _tfs, _nchain) { }

    protected:
        Eigen::ArrayXXd mu;         // extented mass
        Eigen::ArrayXXd eta;        // extented position
        Eigen::ArrayXXd theta;      // extented momentum
        Eigen::ArrayXXd Gamma;      // thermostat force

        void initialize() override;
        void calc_thermo_force(const int& j);
        void calc_cons_quant();
    };

    /*** ================================================== ***/
    /*** NHC Procedure (Local, Side)                        ***/
    /*** ================================================== ***/

    class nhc_local_side : public nhc_local_base {
    public:
        nhc_local_side() = default;
        nhc_local_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            nhc_local_base(_fn_no_ex, _bsp, _sys, _tfs, _nchain)
        {
            labels.push_back("NHC");
            labels.push_back("Local");
            labels.push_back("Side");
        }

    private:
        void implement_one_step() override;
    };

    /*** ================================================== ***/
    /*** NHC Procedure (Local, Middle)                        ***/
    /*** ================================================== ***/

    class nhc_local_middle : public nhc_local_base {
    public:
        nhc_local_middle() = default;
        nhc_local_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const thermo_factor_scheme& _tfs, const int _nchain) :
            nhc_local_base(_fn_no_ex, _bsp, _sys, _tfs, _nchain)
        {
            labels.push_back("NHC");
            labels.push_back("Local");
            labels.push_back("Middle");
        }

    private:
        void implement_one_step() override;
    };

} // !nhc
} // !thermostat
} // !uovie
#endif // !NOSE_HOOVER_CHAIN_H_