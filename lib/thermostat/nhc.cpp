/* Nose-Hoover Chain */
// standard C++ headers
#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>

// uovie headers
#include "thermostat/nhc.h"
#include "model.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** NHC Procedure (Base) Class Member Functions           ***/
    /*** ===================================================== ***/

    void nhc_base::print_ther_proce_title(std::ofstream& chk, std::ofstream& out)
    {
        std::cout << labels[0] + " procedure (" + labels[1] + ", " + labels[2] + ") is running." << std::endl;
        chk << labels[0] + " procedure (" + labels[1] + ", " + labels[2] + "):\n   Time" << "            " << "position"
            << "            " << "momentum" << "            " << "con_ene" << std::endl;
        out << labels[0] + " procedure (" + labels[1] + ", " + labels[2] + "):\n   Time" << "            " << "kin_ene"
            << "            " << "pot_ene" << std::endl;
    }

    void nhc_base::print_ther_proce_data(std::ofstream& chk, std::ofstream& out, double& t)
    {
        chk << std::scientific << std::setprecision(8) << std::setw(20) << t
            << std::setw(20) << q(0) << std::setw(20) << p(0) << std::setw(20) << con_ene << std::endl;
        out << std::scientific << std::setprecision(8) << std::setw(20) << t
            << std::setw(20) << kin_ene << std::setw(20) << pot_ene << std::endl;
    }

    void nhc_base::print_conclusion_info(std::ofstream& chk, std::ofstream& out,
        const std::chrono::duration<double>& time_elap) {
        chk << "\nElapsed time of " + labels[0] + " procedure (" + labels[1] + ", " + labels[2] + "): " << time_elap.count() << std::endl;
        out << "\nElapsed time of " + labels[0] + " procedure (" + labels[1] + ", " + labels[2] + "): " << time_elap.count() << std::endl;
        std::cout << "\nElapsed time of " + labels[0] + " procedure (" + labels[1] + ", " + labels[2] + "): " << time_elap.count() << std::endl;
        std::cout << "\nNormal termination. Congratulations!" << std::endl;
    }

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** NHC Procedure (Global, Base) Class Member Functions   ***/
    /*** ===================================================== ***/

    // initialization
    void nhc_global_base::initialize()
    {
        // resize arrays
        m.resize(dof);
        q.resize(dof);
        p.resize(dof);
        F.resize(dof);

        mu.resize(M);
        eta.resize(M);
        theta.resize(M);
        Gamma.resize(M);

        // initialize positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi) = sys.molecules[mi].atoms[ai].m;
                    q(vi) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }

        // initialize momenta
        std::mt19937 p_mte(36);
        for (auto ri = 0; ri < p.rows(); ri++) {
            std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri)) };
            p(ri) = ndrm(p_mte);
        }

        // initialize extented masses and extented momenta
        std::mt19937 theta_mte(27);
        //if (sys.model_type == "HO") { // simple harmonic forces
        //    model::harmonic_oscilator HO(sys.model_para[0]);
        //    double tau = 1 / HO.ome();
        double tau = 1;
        mu = k * T * pow(tau, 2) * Eigen::ArrayXd::Constant(M, 1);
        for (auto j = 0; j < M; j++) {
            std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu(j)) };
            theta(j) = ndrm(theta_mte);
        }
        //}

        // initialize extented positions
        eta.setZero();
    }

    void nhc_global_base::calc_thermo_force(const int& j)
    {
        if (j == 0)
            Gamma(0) = (p.pow(2) / m).sum() - d * N * k * T;
        else
            Gamma(j) = pow(theta(j - 1), 2) / mu(j - 1) - k * T;
    }

    void nhc_global_base::calc_cons_quant()
    {
        kin_ene = (p.pow(2) / (2 * m)).sum();

        pot_ene = 0;
        if (sys.model_type == "HO") { // harmonic oscillator
            model::harmonic_oscilator HO(sys.model_para[0]);
            pot_ene = 0.5 * pow(HO.ome(), 2) * (m * q.pow(2)).sum();
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ni = 0; ni < N - 1; ni++) {
                for (auto nj = ni + 1; nj < N; nj++) {
                    if (ni == nj) continue;
                    pot_ene += LJ.V(q.block(d * ni, 0, d, 1), q.block(d * nj, 0, d, 1));
                }
            }
        }

        the_ene = pow(theta(0), 2) / (2 * mu(0)) + d * N * k * T * eta(0);
        for (int j = 1; j < M; j++)
            the_ene += pow(theta(j), 2) / (2 * mu(j)) + k * T * eta(j);

        con_ene = kin_ene + pot_ene + the_ene;
    }

    /*** ===================================================== ***/
    /*** NHC Procedure (Global, Side) Class Member Functions   ***/
    /*** ===================================================== ***/

    void nhc_global_side::implement_one_step()
    {
        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta(M - 1) += tmp_delta * Gamma(M - 1) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                    calc_thermo_force(j);
                    theta(j) += tmp_delta * Gamma(j) / 4;
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                }
                eta += tmp_delta * theta / (2 * mu);

                double momen_scale = exp(-1 * tmp_delta * theta(0) / (2 * mu(0)));
                p *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                    calc_thermo_force(j);
                    theta(j) += tmp_delta * Gamma(j) / 4;
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                }
                calc_thermo_force(M - 1);
                theta(M - 1) += tmp_delta * Gamma(M - 1) / 4;
            }
        }

        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / m;
        calc_physic_force();
        p += F * Dt / 2;

        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta(M - 1) += tmp_delta * Gamma(M - 1) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                    calc_thermo_force(j);
                    theta(j) += tmp_delta * Gamma(j) / 4;
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                }
                eta += tmp_delta * theta / (2 * mu);

                double momen_scale = exp(-1 * tmp_delta * theta(0) / (2 * mu(0)));
                p *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                    calc_thermo_force(j);
                    theta(j) += tmp_delta * Gamma(j) / 4;
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (8 * mu(j + 1)));
                }
                calc_thermo_force(M - 1);
                theta(M - 1) += tmp_delta * Gamma(M - 1) / 4;
            }
        }
    }

    /*** ===================================================== ***/
    /*** NHC Procedure (Global, Middle) Class Member Functions ***/
    /*** ===================================================== ***/

    void nhc_global_middle::implement_one_step()
    {
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / (2 * m);

        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta(M - 1) += tmp_delta * Gamma(M - 1) / 2;

                for (int j = M - 2; j >= 0; j--) {
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (4 * mu(j + 1)));
                    calc_thermo_force(j);
                    theta(j) += tmp_delta * Gamma(j) / 2;
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (4 * mu(j + 1)));
                }
                eta += tmp_delta * theta / mu;

                double momen_scale = exp(-1 * tmp_delta * theta(0) / mu(0));
                p *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (4 * mu(j + 1)));
                    calc_thermo_force(j);
                    theta(j) += tmp_delta * Gamma(j) / 2;
                    theta(j) *= exp(-1 * tmp_delta * theta(j + 1) / (4 * mu(j + 1)));
                }
                calc_thermo_force(M - 1);
                theta(M - 1) += tmp_delta * Gamma(M - 1) / 2;
            }
        }

        q += p * Dt / (2 * m);
        calc_physic_force();
        p += F * Dt / 2;
    }

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** NHC Procedure (Local, Base) Class Member Functions    ***/
    /*** ===================================================== ***/

    // initialization
    void nhc_local_base::initialize()
    {
        // resize arrays
        m.resize(dof);
        q.resize(dof);
        p.resize(dof);
        F.resize(dof);

        mu.resize(dof, M);
        eta.resize(dof, M);
        theta.resize(dof, M);
        Gamma.resize(dof, M);

        // initialize positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi) = sys.molecules[mi].atoms[ai].m;
                    q(vi) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }

        // initialize momenta
        std::mt19937 p_mte(36);
        for (auto ri = 0; ri < p.rows(); ri++) {
            std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri)) };
            p(ri) = ndrm(p_mte);
        }

        // initialize extented masses and extented momenta
        std::mt19937 theta_mte(27);
        //if (sys.model_type == "HO") { // simple harmonic forces
        //    model::harmonic_oscilator HO(sys.model_para[0]);
        //   double tau = 1 / HO.ome();
        double tau = 1;
        mu = k * T * pow(tau, 2) * Eigen::ArrayXXd::Constant(d * N, M, 1);
        for (auto ri = 0; ri < theta.rows(); ri++) {
            for (auto ci = 0; ci < theta.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu(ri, ci)) };
                theta(ri, ci) = ndrm(theta_mte);
            }
        }
        //}

        // initialize extented positions
        eta.setZero();

    }

    void nhc_local_base::calc_thermo_force(const int& j)
    {
        if (j == 0)
            Gamma.col(0) = p.pow(2) / m - k * T;
        else
            Gamma.col(j) = theta.col(j - 1).pow(2) / mu.col(j - 1) - k * T;
    }

    void nhc_local_base::calc_cons_quant()
    {
        kin_ene = (p.pow(2) / (2 * m)).sum();

        pot_ene = 0;
        if (sys.model_type == "HO") { // harmonic oscillator
            model::harmonic_oscilator HO(sys.model_para[0]);
            pot_ene = 0.5 * pow(HO.ome(), 2) * (m * q.pow(2)).sum();
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ni = 0; ni < N - 1; ni++) {
                for (auto nj = ni + 1; nj < N; nj++) {
                    if (ni == nj) continue;
                    pot_ene += LJ.V(q.block(d * ni, 0, d, 1), q.block(d * nj, 0, d, 1));
                }
            }
        }

        the_ene = (theta.pow(2) / (2 * mu) + k * T * eta).sum();

        con_ene = kin_ene + pot_ene + the_ene;
    }

    /*** ===================================================== ***/
    /*** NHC Procedure (Local, Side) Class Member Functions    ***/
    /*** ===================================================== ***/

    void nhc_local_side::implement_one_step()
    {
        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.col(M - 1) += tmp_delta * Gamma.col(M - 1) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                    calc_thermo_force(j);
                    theta.col(j) += tmp_delta * Gamma.col(j) / 4;
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                }
                eta += tmp_delta * theta / (2 * mu);

                Eigen::ArrayXd momen_scale = (-1 * tmp_delta * theta.col(0) / (2 * mu.col(0))).exp();
                p *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                    calc_thermo_force(j);
                    theta.col(j) += tmp_delta * Gamma.col(j) / 4;
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                }
                calc_thermo_force(M - 1);
                theta.col(M - 1) += tmp_delta * Gamma.col(M - 1) / 4;
            }
        }

        calc_physic_force();
        p += Dt * F / 2;
        q += Dt * p / m;
        calc_physic_force();
        p += Dt * F / 2;

        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.col(M - 1) += tmp_delta * Gamma.col(M - 1) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                    calc_thermo_force(j);
                    theta.col(j) += tmp_delta * Gamma.col(j) / 4;
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                }
                eta += tmp_delta * theta / (2 * mu);

                Eigen::ArrayXd momen_scale = (-1 * tmp_delta * theta.col(0) / (2 * mu.col(0))).exp();
                p *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                    calc_thermo_force(j);
                    theta.col(j) += tmp_delta * Gamma.col(j) / 4;
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (8 * mu.col(j + 1))).exp();
                }
                calc_thermo_force(M - 1);
                theta.col(M - 1) += tmp_delta * Gamma.col(M - 1) / 4;
            }
        }
    }

    /*** ===================================================== ***/
    /*** NHC Procedure (Local, Middle) Class Member Functions  ***/
    /*** ===================================================== ***/

    void nhc_local_middle::implement_one_step()
    {
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / (2 * m);

        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.col(M - 1) += tmp_delta * Gamma.col(M - 1) / 2;

                for (int j = M - 2; j >= 0; j--) {
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (4 * mu.col(j + 1))).exp();
                    calc_thermo_force(j);
                    theta.col(j) += tmp_delta * Gamma.col(j) / 2;
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (4 * mu.col(j + 1))).exp();
                }
                eta += tmp_delta * theta / mu;

                Eigen::ArrayXd momen_scale = (-1 * tmp_delta * theta.col(0) / mu.col(0)).exp();
                p *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (4 * mu.col(j + 1))).exp();
                    calc_thermo_force(j);
                    theta.col(j) += tmp_delta * Gamma.col(j) / 2;
                    theta.col(j) *= (-1 * tmp_delta * theta.col(j + 1) / (4 * mu.col(j + 1))).exp();
                }
                calc_thermo_force(M - 1);
                theta.col(M - 1) += tmp_delta * Gamma.col(M - 1) / 2;
            }
        }

        q += p * Dt / (2 * m);
        calc_physic_force();
        p += F * Dt / 2;
    }

} // !nhc
} // !thermostat
} // !uovie