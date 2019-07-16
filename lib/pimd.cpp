// standard C++ headers
#include <iostream>
#include <vector>
#include <iomanip>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "phy_const.h"
#include "pimd.h"
#include "model.h"

namespace uovie {
namespace pimd {

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** PIMD Procedure (Base) Class Member Functions          ***/
    /*** ===================================================== ***/

    // staging transformation
    void pimd_base::stag_trans()
    {
        r.col(0) = q.col(0);
        for (int ci = 1; ci < r.cols() - 1; ci++)
            r.col(ci) = -(1 / (ci + 1.0)) * q.col(0) + q.col(ci) - (ci / (ci + 1.0)) * q.col(ci + 1);
        r.col(r.cols() - 1) = q.col(r.cols() - 1) - q.col(0);
    }

    // inverse staging transformation
    void pimd_base::inve_stag_trans()
    {
        q.col(q.cols() - 1) = r.col(q.cols() - 1) + r.col(0);
        for (int ci = q.cols() - 2; ci > 0; ci--)
            q.col(ci) = (1 / (ci + 1.0)) * r.col(0) + r.col(ci) + (ci / (ci + 1.0)) * q.col(ci + 1);
        q.col(0) = r.col(0);
    }

    // calculate physical forces
    void pimd_base::calc_physic_force()
    {
        Eigen::ArrayXXd pV_pq = Eigen::ArrayXXd::Zero(dof, nbead); // \frac{\partial V}{\partial q}
        Eigen::ArrayXXd pV_pr = Eigen::ArrayXXd::Zero(dof, nbead); // \frac{\partial V}{\partial r}

        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            pV_pq = m * pow(HO.ome(), 2) * q;
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ci = 0; ci < pV_pq.cols(); ci++) {
                for (auto ni = 0; ni < N; ni++) {
                    for (auto nj = 0; nj < N; nj++) {
                        if (ni == nj) continue;
                        pV_pq.block(d * ni, ci, d, 1) -= LJ.F(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
                    }
                }
            }
        }

        for (auto cj = 0; cj < pV_pr.cols(); cj++)
            pV_pr.col(0) += pV_pq.col(cj);
        for (auto ci = 1; ci < pV_pr.cols(); ci++)
            pV_pr.col(ci) = pV_pq.col(ci) + ((ci - 1.0) / ci) * pV_pr.col(ci - 1);
        F = -1 * m_bar * pow(fic_omega, 2) * r - pV_pr / nbead;

    }

    void pimd_base::calc_prim_estor()
    {
        prim_kin_estor = 0;
        for (int bi = 0; bi < nbead - 1; bi++)
            prim_kin_estor += (m.col(bi) * (q.col(bi + 1) - q.col(bi)).pow(2)).sum();
        prim_kin_estor += (m.col(nbead - 1) * (q.col(0) - q.col(nbead - 1)).pow(2)).sum();
        prim_kin_estor *= -1 * nbead / (2 * pow(h_bar * beta, 2));
        prim_kin_estor += dof * nbead / (2 * beta);

        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            pot_ene = 0.5 * pow(HO.ome(), 2) * (m * q.pow(2)).sum();
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            pot_ene = 0;
            for (auto ci = 0; ci < nbead; ci++)
                for (auto ni = 0; ni < N - 1; ni++)
                    for (auto nj = ni + 1; nj < N; nj++)
                        pot_ene += LJ.V(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
        }
        prim_pot_estor = pot_ene / nbead;

        //prim_pres_estor = N*nbead/(beta * V)
    }

    void pimd_base::implement()
    {
        initialize();

        double t = 0;
        int ctr = 0;

        std::ofstream chk, out;
        chk.open(fn_no_ex + ".chk");
        out.open(fn_no_ex + ".out");

        print_pimd_proce_title(chk, out);
        calc_prim_estor();
        print_pimd_proce_data(chk, out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            implement_one_step();
            if (ctr == bsp.data_coll_peri) {
                calc_prim_estor();
                print_pimd_proce_data(chk, out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        print_conclusion_info(chk, out, time_elapsed);
        chk.close();
        out.close();
    }

    void pimd_base::print_pimd_proce_title(std::ofstream& chk, std::ofstream& out)
    {
        std::cout << "\nPIMD procedure via " + labels[0] + " (" + labels[1] + ") is running." << std::endl;
        chk << "# PIMD procedure via " + labels[0] + " (" + labels[1] + "):\n" << "#          Time" << "              position"
            << "            momentum" << std::endl;
        out << "# PIMD procedure via " + labels[0] + " (" + labels[1] + "):\n" << "#          Time" << "           prim_kin_estor"
            << "      prim_pot_estor" << std::endl;
    }

    void pimd_base::print_pimd_proce_data(std::ofstream& chk, std::ofstream& out, double& t)
    {
        chk << std::scientific << std::setprecision(8) << std::setw(20) << t * uovie::phy_const::a_u_time * 1e15
            << std::setw(20) << q(0, 0) << std::setw(20) << s(0, 0) << std::endl;
        out << std::scientific << std::setprecision(8) << std::setw(20) << t * uovie::phy_const::a_u_time * 1e15
            << std::setw(20) << prim_kin_estor << std::setw(20) << prim_pot_estor << std::endl;
    }

    void pimd_base::print_conclusion_info(std::ofstream& chk, std::ofstream& out,
        const std::chrono::duration<double>& time_elap) {
        chk << "\n# Elapsed time of PIMD procedure via " + labels[0] + " (" + labels[1] + "): " << time_elap.count() << std::endl;
        out << "\n# Elapsed time of PIMD procedure via " + labels[0] + " (" + labels[1] + "): " << time_elap.count() << std::endl;
        std::cout << "\nElapsed time of PIMD procedure via " + labels[0] + " (" + labels[1] + "): " << time_elap.count() << std::endl;
        std::cout << "\nNormal termination. Congratulations!" << std::endl;
    }

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** PIMD Procedure via LD (Base) Class Member Functions   ***/
    /*** ===================================================== ***/

    // initialization
    void pimd_via_ld_base::initialize()
    {
        // resize arrays
        m.resize(dof, nbead);
        q.resize(dof, nbead);
        m_tilde.resize(dof, nbead);
        r.resize(dof, nbead);
        s.resize(dof, nbead);
        F.resize(dof, nbead);
        nrand.resize(dof, nbead);

        // initialize cartisian positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi, 0) = sys.molecules[mi].atoms[ai].m;
                    q(vi, 0) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }
        for (auto ci = 1; ci < q.cols(); ci++) {
            q.col(ci) = q.col(0);
            m.col(ci) = m.col(0);
        }

        // initialize staging transformed masses
        m_tilde.col(0) = m.col(0);
        for (auto ci = 1; ci < m_tilde.cols(); ci++)
            m_tilde.col(ci) = ((ci + 1.0) / ci) * m.col(ci);

        m_bar = m_tilde;
        m_bar.col(0).setZero();

        // initialize staging transformed positions
        stag_trans();

        // initialize fictition momenta
        std::random_device rd;
        std::mt19937 s_mte(rd());
        for (auto ri = 0; ri < s.rows(); ri++) {
            for (auto ci = 0; ci < s.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri, ci)) };
                s(ri, ci) = ndrm(s_mte);
            }
        }

        // set random-number generators
        std::mt19937 tmp_mte;
        std::normal_distribution<double> tmp_nd{ 0, 1 };
        for (int i = 0; i < dof * nbead; i++) {
            mtes.push_back(tmp_mte);
            mtes[i].seed(rd());
            nds.push_back(tmp_nd);
        }

    }

    void pimd_via_ld_base::upd_rand_num()
    {
        for (int ci = 0; ci < nbead; ci++)
            for (int ri = 0; ri < dof; ri++)
                nrand(ri, ci) = nds[ri + ci * dof](mtes[ri + ci * dof]);
    }

    /*** ===================================================== ***/
    /*** PIMD Procedure via LD (Side) Class Member Functions   ***/
    /*** ===================================================== ***/

    void pimd_via_ld_side::implement_one_step()
    {
        upd_rand_num();
        s = c1 * s + c2 * sqrt(1 / beta) * m_tilde.sqrt() * nrand;
        calc_physic_force();
        s += F * Dt / 2;
        r += s * Dt / m_tilde;
        inve_stag_trans();
        calc_physic_force();
        s += F * Dt / 2;
        upd_rand_num();
        s = c1 * s + c2 * sqrt(1 / beta) * m_tilde.sqrt() * nrand;
    }

    /*** ===================================================== ***/
    /*** PIMD Procedure via LD (Middle) Class Member Functions ***/
    /*** ===================================================== ***/

    void pimd_via_ld_middle::implement_one_step()
    {
        calc_physic_force();
        s += F * Dt / 2;
        r += s * Dt / (2 * m_tilde);

        upd_rand_num();
        s = c1 * s + c2 * sqrt(1 / beta) * m_tilde.sqrt() * nrand;

        r += s * Dt / (2 * m_tilde);
        inve_stag_trans();
        calc_physic_force();
        s += F * Dt / 2;
    }

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** PIMD Procedure via AT (Base) Class Member Functions   ***/
    /*** ===================================================== ***/

    // initialization
    void pimd_via_at_base::initialize()
    {
        // resize arrays
        m.resize(dof, nbead);
        q.resize(dof, nbead);
        m_tilde.resize(dof, nbead);
        r.resize(dof, nbead);
        s.resize(dof, nbead);
        F.resize(dof, nbead);
        nrand.resize(dof, nbead);
        urand.resize(dof, nbead);

        // initialize cartisian positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi, 0) = sys.molecules[mi].atoms[ai].m;
                    q(vi, 0) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }
        for (auto ci = 1; ci < q.cols(); ci++) {
            q.col(ci) = q.col(0);
            m.col(ci) = m.col(0);
        }

        // initialize staging transformed masses
        m_tilde.col(0) = m.col(0);
        for (auto ci = 1; ci < m_tilde.cols(); ci++)
            m_tilde.col(ci) = ((ci + 1.0) / ci) * m.col(ci);

        m_bar = m_tilde;
        m_bar.col(0).setZero();

        // initialize staging transformed positions
        stag_trans();

        // initialize fictition momenta
        std::random_device rd;
        std::mt19937 s_mte(rd());
        for (auto ri = 0; ri < s.rows(); ri++) {
            for (auto ci = 0; ci < s.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri, ci)) };
                s(ri, ci) = ndrm(s_mte);
            }
        }

        // set normal-distributed random-number generators
        std::mt19937 tmp_nmte;
        std::normal_distribution<double> tmp_nd{ 0, 1 };
        for (int i = 0; i < dof * nbead; i++) {
            nmtes.push_back(tmp_nmte);
            nmtes[i].seed(rd());
            nds.push_back(tmp_nd);
        }

        // set uniform-distributed random-number generators
        std::mt19937 tmp_umte;
        std::uniform_real_distribution<double> tmp_ud{ 0, 1 };
        for (int i = 0; i < dof * nbead; i++) {
            umtes.push_back(tmp_umte);
            umtes[i].seed(rd());
            urds.push_back(tmp_ud);
        }
    }

    void pimd_via_at_base::upd_uni_rand_num()
    {
        for (int ci = 0; ci < nbead; ci++)
            for (int ri = 0; ri < dof; ri++)
                urand(ri, ci) = urds[ri + ci * dof](umtes[ri + ci * dof]);
    }

    /*** ===================================================== ***/
    /*** PIMD Procedure via AT (Side) Class Member Functions   ***/
    /*** ===================================================== ***/

    void pimd_via_at_side::implement_one_step()
    {
        upd_uni_rand_num();
        for (int ci = 0; ci < nbead; ci++) {
            for (int ri = 0; ri < dof; ri++) {
                if (urand(ri, ci) < cri) {
                    nrand(ri, ci) = nds[ri + ci * dof](nmtes[ri + ci * dof]);
                    s(ri, ci) = sqrt(1 / beta) * m_tilde(ri, ci) * nrand(ri, ci);
                }
            }
        }

        calc_physic_force();
        s += F * Dt / 2;
        r += s * Dt / m;
        inve_stag_trans();
        calc_physic_force();
        s += F * Dt / 2;

        upd_uni_rand_num();
        for (int ci = 0; ci < nbead; ci++) {
            for (int ri = 0; ri < dof; ri++) {
                if (urand(ri, ci) < cri) {
                    nrand(ri, ci) = nds[ri + ci * dof](nmtes[ri + ci * dof]);
                    s(ri, ci) = sqrt(1 / beta) * m_tilde(ri, ci) * nrand(ri, ci);
                }
            }
        }
    }

    /*** ===================================================== ***/
    /*** PIMD Procedure via AT (Middle) Class Member Functions ***/
    /*** ===================================================== ***/

    void pimd_via_at_middle::implement_one_step()
    {
        calc_physic_force();
        s += F * Dt / 2;
        r += s * Dt / (2 * m);

        upd_uni_rand_num();
        for (int ci = 0; ci < nbead; ci++) {
            for (int ri = 0; ri < dof; ri++) {
                if (urand(ri, ci) < cri) {
                    nrand(ri, ci) = nds[ri + ci * dof](nmtes[ri + ci * dof]);
                    s(ri, ci) = sqrt(1 / beta) * m_tilde(ri, ci) * nrand(ri, ci);
                }
            }
        }

        r += s * Dt / (2 * m);
        inve_stag_trans();
        calc_physic_force();
        s += F * Dt / 2;
    }

    //-----------------------------------------------------------//

    /*** ===================================================== ***/
    /*** PIMD Procedure via NHC (Base) Class Member Functions  ***/
    /*** ===================================================== ***/

    // initialization
    void pimd_via_nhc_base::initialize()
    {
        // resize arrays
        m.resize(dof, nbead);
        q.resize(dof, nbead);
        m_tilde.resize(dof, nbead);
        r.resize(dof, nbead);
        s.resize(dof, nbead);
        F.resize(dof, nbead);

        mu.resize(dof * M, nbead);
        eta.resize(dof * M, nbead);
        theta.resize(dof * M, nbead);
        Gamma.resize(dof * M, nbead);

        // initialize cartisian positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi, 0) = sys.molecules[mi].atoms[ai].m;
                    q(vi, 0) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }
        for (auto ci = 1; ci < q.cols(); ci++) {
            q.col(ci) = q.col(0);
            m.col(ci) = m.col(0);
        }

        // initialize staging transformed masses
        m_tilde.col(0) = m.col(0);
        for (auto ci = 1; ci < m_tilde.cols(); ci++)
            m_tilde.col(ci) = ((ci + 1.0) / ci) * m.col(ci);

        m_bar = m_tilde;
        m_bar.col(0).setZero();

        // initialize staging transformed positions
        stag_trans();

        // initialize fictition momenta
        std::random_device rd;
        std::mt19937 s_mte(rd());
        for (auto ri = 0; ri < s.rows(); ri++) {
            for (auto ci = 0; ci < s.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri, ci)) };
                s(ri, ci) = ndrm(s_mte);
            }
        }

        // initialize extented masses and extented momenta
        std::mt19937 theta_mte(rd());
        double tau = 1;
        if (sys.model_type == "LJ")
            tau = 5.16767166e4;
        mu = k * T * pow(tau, 2) * Eigen::ArrayXXd::Constant(d * N * M, nbead, 1);
        for (auto ri = 0; ri < theta.rows(); ri++) {
            for (auto ci = 0; ci < theta.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu(ri, ci)) };
                theta(ri, ci) = ndrm(theta_mte);
            }
        }

        // initialize extented positions
        eta.setZero();

    }

    void pimd_via_nhc_base::calc_thermo_force(const int& j)
    {
        if (j == 0)
            Gamma.block(0, 0, dof, nbead) = s.pow(2) / m_tilde - k * T;
        else
            Gamma.block(j * dof, 0, dof, nbead) = theta.block((j - 1) * dof, 0, dof, nbead).pow(2)
            / mu.block((j - 1) * dof, 0, dof, nbead) - k * T;
    }

    void pimd_via_nhc_base::calc_cons_quant()
    {
        kin_ene = (s.pow(2) / (2 * m_tilde)).sum();

        if (sys.model_type == "HO") { // harmonic oscillator
            model::harmonic_oscilator HO(sys.model_para[0]);
            pot_ene = 0.5 * pow(HO.ome(), 2) * (m * q.pow(2)).sum();
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            pot_ene = 0;
            for (auto ci = 0; ci < nbead; ci++)
                for (auto ni = 0; ni < N - 1; ni++)
                    for (auto nj = ni + 1; nj < N; nj++)
                        pot_ene += LJ.V(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
        }

        double the_ene = 0;
        Eigen::ArrayXXd the_ene_arr = Eigen::ArrayXXd::Zero(dof, nbead);
        for (int j = 0; j < M; j++)
            the_ene_arr += theta.block(j * dof, 0, dof, nbead).pow(2) / (2 * mu.block(j * dof, 0, dof, nbead))
            + k * T * eta.block(j * dof, 0, dof, nbead);
        the_ene = the_ene_arr.sum();

        con_ene = kin_ene + 0.5 * pow(fic_omega, 2) * (m_bar * r.pow(2)).sum() + pot_ene / nbead + the_ene;
    }

    void pimd_via_nhc_base::implement()
    {
        initialize();

        double t = 0;
        int ctr = 0;

        std::ofstream chk, out;
        chk.open(fn_no_ex + ".chk");
        out.open(fn_no_ex + ".out");

        print_pimd_proce_title(chk, out);
        calc_prim_estor();
        calc_cons_quant();
        print_pimd_proce_data(chk, out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            implement_one_step();
            if (ctr == bsp.data_coll_peri) {
                calc_prim_estor();
                calc_cons_quant();
                print_pimd_proce_data(chk, out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        print_conclusion_info(chk, out, time_elapsed);
        chk.close();
        out.close();
    }

    void pimd_via_nhc_base::print_pimd_proce_title(std::ofstream& chk, std::ofstream& out)
    {
        std::cout << "\nPIMD procedure via " + labels[0] + " (" + labels[1] + ") is running." << std::endl;
        chk << "# PIMD procedure via " + labels[0] + " (" + labels[1] + "):\n" << "#          Time" << "              position"
            << "            momentum" << "            con_ene" << std::endl;
        out << "# PIMD procedure via " + labels[0] + " (" + labels[1] + "):\n" << "#          Time" << "           prim_kin_estor"
            << "      prim_pot_estor" << std::endl;
    }

    void pimd_via_nhc_base::print_pimd_proce_data(std::ofstream& chk, std::ofstream& out, double& t)
    {
        chk << std::scientific << std::setprecision(8) << std::setw(20) << t * uovie::phy_const::a_u_time * 1e15
            << std::setw(20) << q(0, 0) << std::setw(20) << s(0, 0) << std::setw(20) << con_ene << std::endl;
        out << std::scientific << std::setprecision(8) << std::setw(20) << t * uovie::phy_const::a_u_time * 1e15
            << std::setw(20) << prim_kin_estor << std::setw(20) << prim_pot_estor << std::endl;
    }

    /*** ===================================================== ***/
    /*** PIMD Procedure via NHC (Side) Class Member Functions  ***/
    /*** ===================================================== ***/

    void pimd_via_nhc_side::implement_one_step()
    {
        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 4;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                eta += tmp_delta * theta / (2 * mu);

                Eigen::ArrayXXd momen_scale = (-1 * tmp_delta * theta.block(0, 0, dof, nbead)
                    / (2 * mu.block(0, 0, dof, nbead))).exp();
                s *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 4;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 4;
            }
        }

        calc_physic_force();
        s += Dt * F / 2;
        r += Dt * s / m_tilde;
        inve_stag_trans();
        calc_physic_force();
        s += Dt * F / 2;

        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 4;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                eta += tmp_delta * theta / (2 * mu);

                Eigen::ArrayXXd momen_scale = (-1 * tmp_delta * theta.block(0, 0, dof, nbead)
                    / (2 * mu.block(0, 0, dof, nbead))).exp();
                s *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 4;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 4;
            }
        }
    }

    /*** ====================================================== ***/
    /*** PIMD Procedure via NHC (Middle) Class Member Functions ***/
    /*** ====================================================== ***/

    void pimd_via_nhc_middle::implement_one_step()
    {
        calc_physic_force();
        s += Dt * F / 2;
        r += Dt * s / (2 * m_tilde);

        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 2;

                for (int j = M - 2; j >= 0; j--) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (4 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 2;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (4 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                eta += tmp_delta * theta / mu;

                Eigen::ArrayXXd momen_scale = (-1 * tmp_delta * theta.block(0, 0, dof, nbead)
                    / mu.block(0, 0, dof, nbead)).exp();
                s *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (4 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 2;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (4 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 2;
            }
        }

        r += Dt * s / (2 * m_tilde);
        inve_stag_trans();
        calc_physic_force();
        s += Dt * F / 2;
    }

} // !pimd
} // !uovie