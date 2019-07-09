/* Andersen Thermostat */

// standard C++ headers
#include <iostream>
#include <iomanip>
#include <chrono>

// uovie headers
#include "thermostat/at.h"
#include "model.h"

namespace uovie {
namespace thermostat {
namespace at {

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** AT Procedure (Side) Class Member Functions         ***/
    /*** ================================================== ***/

    // initialization
    void at_procedure_side::initialize()
    {
        // resize arrays
        m.resize(dof);
        q.resize(dof);
        p.resize(dof);
        F.resize(dof);
        nrand.resize(dof);
        urand.resize(dof);

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

        // set normal-distributed random-number generators
        std::random_device rd;
        std::mt19937 tmp_nmte;
        std::normal_distribution<double> tmp_nd{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            nmtes.push_back(tmp_nmte);
            nmtes[i].seed(rd());
            nds.push_back(tmp_nd);
        }

        // set uniform-distributed random-number generators
        std::mt19937 tmp_umte;
        std::uniform_real_distribution<double> tmp_ud{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            umtes.push_back(tmp_umte);
            umtes[i].seed(rd());
            urds.push_back(tmp_ud);
        }
        
    }

    // calculate physical forces
    void at_procedure_side::calc_physic_force()
    {
        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            F = -1 * m * pow(HO.ome(), 2) * q;
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            F.setZero();
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ni = 0; ni < N; ni++) {
                for (auto nj = 0; nj < N; nj++) {
                    if (ni == nj) continue;
                    F.block(d * ni, 0, d, 1) += LJ.F(q.block(d * ni, 0, d, 1), q.block(d * nj, 0, d, 1));
                }
            }
        }
    }

    void at_procedure_side::upd_uni_rand_num()
    {
        for (int i = 0; i < dof; i++)
            urand(i) = urds[i](umtes[i]);
    }

    void at_procedure_side::print_at_procedure_title(std::ofstream& out) {
        std::cout << "\nAT Procedure (Side):\n   Time" << "            " << "position"
            << "            " << "momentum";
        out << "\nAT Procedure (Side):\n   Time" << "            " << "position"
            << "            " << "momentum";
    }

    void at_procedure_side::print_at_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::scientific << std::setw(20) << q(0) << std::setw(20) << p(0);
        out << std::scientific << std::setw(20) << q(0) << std::setw(20) << p(0);
    }

    void at_procedure_side::implement_one_step()
    {
        upd_uni_rand_num();
        for (int i = 0; i < dof; i++) {
            if (urand(i) < cri) {
                nrand(i) = nds[i](nmtes[i]);
                p(i) = sqrt(1 / beta) * m(i) * nrand(i);
            }
        }
        
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / m;
        calc_physic_force();
        p += F * Dt / 2;

        upd_uni_rand_num();
        for (int i = 0; i < dof; i++) {
            if (urand(i) < cri) {
                nrand(i) = nds[i](nmtes[i]);
                p(i) = sqrt(1 / beta) * m(i) * nrand(i);
            }
        }
    }

    void at_procedure_side::implement() {
        for (double t = 0; t <= bsp.run_time; t += Dt) {
            implement_one_step();
            t += Dt;
        }
    }

    void at_procedure_side::implement(std::ofstream& out)
    {
        initialize();

        double t = 0;
        int ctr = 0;

        print_at_procedure_title(out);
        print_at_procedure_data(out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            // NHC numerical evolution
            implement_one_step();

            if (ctr == bsp.data_coll_peri) {
                print_at_procedure_data(out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of AT procedure (s): " << time_elapsed.count() << std::endl;
    }

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** AT Procedure (Middle) Class Member Functions       ***/
    /*** ================================================== ***/

    // initialization
    void at_procedure_middle::initialize()
    {
        // resize arrays
        m.resize(dof);
        q.resize(dof);
        p.resize(dof);
        F.resize(dof);
        nrand.resize(dof);
        urand.resize(dof);

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

        // set normal-distributed random-number generators
        std::random_device rd;
        std::mt19937 tmp_nmte;
        std::normal_distribution<double> tmp_nd{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            nmtes.push_back(tmp_nmte);
            nmtes[i].seed(rd());
            nds.push_back(tmp_nd);
        }

        // set uniform-distributed random-number generators
        std::mt19937 tmp_umte;
        std::uniform_real_distribution<double> tmp_ud{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            umtes.push_back(tmp_umte);
            umtes[i].seed(rd());
            urds.push_back(tmp_ud);
        }

    }

    // calculate physical forces
    void at_procedure_middle::calc_physic_force()
    {
        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            F = -1 * m * pow(HO.ome(), 2) * q;
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            F.setZero();
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ni = 0; ni < N; ni++) {
                for (auto nj = 0; nj < N; nj++) {
                    if (ni == nj) continue;
                    F.block(d * ni, 0, d, 1) += LJ.F(q.block(d * ni, 0, d, 1), q.block(d * nj, 0, d, 1));
                }
            }
        }
    }

    void at_procedure_middle::upd_uni_rand_num()
    {
        for (int i = 0; i < dof; i++)
            urand(i) = urds[i](umtes[i]);
    }

    void at_procedure_middle::print_at_procedure_title(std::ofstream& out) {
        std::cout << "\nAT Procedure (Side):\n   Time" << "            " << "position"
            << "            " << "momentum";
        out << "\nAT Procedure (Side):\n   Time" << "            " << "position"
            << "            " << "momentum";
    }

    void at_procedure_middle::print_at_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::scientific << std::setw(20) << q(0) << std::setw(20) << p(0);
        out << std::scientific << std::setw(20) << q(0) << std::setw(20) << p(0);
    }

    void at_procedure_middle::implement_one_step()
    {
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / (2 * m);

        upd_uni_rand_num();
        for (int i = 0; i < dof; i++) {
            if (urand(i) < cri) {
                nrand(i) = nds[i](nmtes[i]);
                p(i) = sqrt(1 / beta) * m(i) * nrand(i);
            }
        }

        q += p * Dt / (2 * m);
        calc_physic_force();
        p += F * Dt / 2;
    }

    void at_procedure_middle::implement() {
        for (double t = 0; t <= bsp.run_time; t += Dt) {
            implement_one_step();
            t += Dt;
        }
    }

    void at_procedure_middle::implement(std::ofstream& out)
    {
        initialize();

        double t = 0;
        int ctr = 0;

        print_at_procedure_title(out);
        print_at_procedure_data(out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            // NHC numerical evolution
            implement_one_step();

            if (ctr == bsp.data_coll_peri) {
                print_at_procedure_data(out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of AT procedure (s): " << time_elapsed.count() << std::endl;
    }

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** AT Procedure for PIMD Class Member Functions      ***/
    /*** ================================================== ***/

    // initialization
    void at_procedure_for_pimd::initialize()
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

        kine_energy.resize(dof, nbead);
        pote_energy.resize(dof, nbead);

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
        std::mt19937 s_mte(36);
        for (auto ri = 0; ri < s.rows(); ri++) {
            for (auto ci = 0; ci < s.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri, ci)) };
                s(ri, ci) = ndrm(s_mte);
            }
        }

        // set normal-distributed random-number generators
        std::random_device rd;
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

    // staging transformation
    void at_procedure_for_pimd::stag_trans()
    {
        r.col(0) = q.col(0);
        for (int ci = 1; ci < r.cols() - 1; ci++)
            r.col(ci) = -(1 / (ci + 1.0)) * q.col(0) + q.col(ci) - (ci / (ci + 1.0)) * q.col(ci + 1);
        r.col(r.cols() - 1) = q.col(r.cols() - 1) - q.col(0);
    }

    // inverse staging transformation
    void at_procedure_for_pimd::inve_stag_trans()
    {
        q.col(q.cols() - 1) = r.col(q.cols() - 1) + r.col(0);
        for (int ci = q.cols() - 2; ci > 0; ci--)
            q.col(ci) = (1 / (ci + 1.0)) * r.col(0) + r.col(ci) + (ci / (ci + 1.0)) * q.col(ci + 1);
        q.col(0) = r.col(0);
    }

    // calculate physical forces
    void at_procedure_for_pimd::calc_physic_force()
    {
        Eigen::ArrayXXd pV_pq = Eigen::ArrayXXd::Zero(d * N, nbead); // \frac{\partial V}{\partial q}
        Eigen::ArrayXXd pV_pr = Eigen::ArrayXXd::Zero(d * N, nbead); // \frac{\partial V}{\partial r}

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

    void at_procedure_for_pimd::upd_uni_rand_num()
    {
        for (int ri = 0; ri < dof; ri++)
            for (int ci = 0; ci < nbead; ci++)
                urand(ri, ci) = urds[ri * dof + ci](umtes[ri * dof + ci]);
    }

    void at_procedure_for_pimd::calc_prim_estor()
    {
        prim_kine_estor = 0;
        for (int bi = 0; bi < nbead - 1; bi++)
            prim_kine_estor += (m.col(bi) * (q.col(bi + 1) - q.col(bi)).pow(2)).sum();
        prim_kine_estor += (m.col(nbead - 1) * (q.col(0) - q.col(nbead - 1)).pow(2)).sum();
        prim_kine_estor *= -1 * nbead / (2 * pow(h_bar * beta, 2));
        prim_kine_estor += dof * nbead / (2 * beta);

        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            pote_energy = 0.5 * m * pow(HO.ome(), 2) * q.pow(2);
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ci = 0; ci < pote_energy.cols(); ci++) {
                for (auto ni = 0; ni < N; ni++) {
                    for (auto nj = 0; nj < N; nj++) {
                        if (ni == nj) continue;
                        pote_energy.block(d * ni, ci, d, 1) += LJ.V(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
                    }
                }
            }
            pote_energy /= 2;
        }
        prim_pote_estor = pote_energy.sum() / nbead;

        //prim_pres_estor = N*nbead/(beta * V)
    }

    void at_procedure_for_pimd::print_at_procedure_title(std::ofstream& out) {
        std::cout << "\nAT Procedure for PIMD:\n   Time" << "            " << "position" << "            " << "momentum"
            << "         " << "prim_kine_estor" << "     " << "prim_pote_estor";
        out << "\nAT Procedure for PIMD:\n   Time" << "            " << "position" << "            " << "momentum"
            << "         " << "prim_kine_estor" << "     " << "prim_pote_estor";
    }

    void at_procedure_for_pimd::print_at_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::scientific << std::setw(20) << q(0, 0) << std::setw(20) << s(0, 0) << std::setw(20)
            << prim_kine_estor << std::setw(20) << prim_pote_estor;
        out << std::scientific << std::setw(20) << q(0, 0) << std::setw(20) << s(0, 0) << std::setw(20)
            << prim_kine_estor << std::setw(20) << prim_pote_estor;
    }

    void at_procedure_for_pimd::implement_one_step() {
        calc_physic_force();
        s += F * Dt / 2;
        r += s * Dt / (2 * m);

        upd_uni_rand_num();
        for (int ri = 0; ri < dof; ri++) {
            for (int ci = 0; ci < nbead; ci++) {
                if (urand(ri, ci) < cri) {
                    nrand(ri, ci) = nds[ri * dof + ci](nmtes[ri * dof + ci]);
                    s(ri, ci) = sqrt(1 / beta) * m_tilde(ri, ci) * nrand(ri, ci);
                }
            }
        }

        r += s * Dt / (2 * m);
        inve_stag_trans();
        calc_physic_force();
        s += F * Dt / 2;
    }

    void at_procedure_for_pimd::implement() {
        for (double t = 0; t <= bsp.run_time; t += Dt) {
            implement_one_step();
            t += Dt;
        }
    }

    void at_procedure_for_pimd::implement(std::ofstream& out)
    {
        initialize();

        double t = 0;
        int ctr = 0;

        print_at_procedure_title(out);
        calc_prim_estor();
        print_at_procedure_data(out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            // NHC numerical evolution
            implement_one_step();

            if (ctr == bsp.data_coll_peri) {
                calc_prim_estor();
                print_at_procedure_data(out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of AT procedure for PIMD (s): " << time_elapsed.count() << std::endl;
    }

} // !ld
} // !thermostat
} // !uovie