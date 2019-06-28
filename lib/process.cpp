#include <iostream>
#include <sstream>
#include <cassert>
#include "../include/process.h"
#include "../include/simu_para.h"
#include "../include/phy_const.h"
#include "../include/atom_data.h"

using namespace uovie::Global;

void process::open(const std::string& file_name) {
    // check the extension
    filename = file_name;
    if (filename.rfind(".vie") == std::string::npos)
        throw std::invalid_argument("Only .vie file is a valid input file.");

    in.open(filename);

    if (in.fail()) {
        std::cout << "Can not open the file " << filename << '.' << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string filename_no_extension;
    for (auto it = filename.begin(); it != filename.end() - 4; it++)
        filename_no_extension += *it;

    // open output file
    out.open(filename_no_extension + ".out");
}

void process::read() {

    constexpr double angstrom_to_bohr = 1e-10 / phy_const::a_u_length;

    std::cout << "Read Infomation from " << filename << std::endl;
    in >> job >> bsp.run_time >> bsp.time_step_size >> bsp.data_coll_peri;
    in >> sys.dimension >> sys.volume >> sys.temperature >> sys.pressure;

    std::string line;
    
    int mol_ctr = 0;
    molecule tmp_mole;
    
    while (getline(in, line)) {
        std::string str;
        if (line.size() >= 5) {
            for (int a = 0; a < 5; a++)
                str.push_back(line[a]);
        }
        if (str == "*****") {
            int num_atom;
            in >> num_atom;
            sys.molecules.push_back(tmp_mole);
            for (int i = 0; i < num_atom; i++) {
                atom tmp_atom;
                sys.molecules[mol_ctr].atoms.push_back(tmp_atom);
                in >> sys.molecules[mol_ctr].atoms[i].symbol;
                for (int d = 0; d < sys.dimension; d++) {
                    double tmp_q, tmp_p, tmp_F;
                    sys.molecules[mol_ctr].atoms[i].q.push_back(tmp_q);
                    sys.molecules[mol_ctr].atoms[i].p.push_back(tmp_p);
                    sys.molecules[mol_ctr].atoms[i].F.push_back(tmp_F);
                    in >> sys.molecules[mol_ctr].atoms[i].q[d];
                    sys.molecules[mol_ctr].atoms[i].p[d] = 1;
                    sys.molecules[mol_ctr].atoms[i].F[d] = 0;
                }
            }
            // identify atomic numbers and assign atomic masses
            for (auto atom_it = sys.molecules[mol_ctr].atoms.begin(); atom_it != sys.molecules[mol_ctr].atoms.end(); atom_it++) {
                for (auto ele_it = atom_data::element.begin(); ele_it != atom_data::element.end(); ele_it++) {
                    if ((*atom_it).symbol == *ele_it) {
                        (*atom_it).atomic_number = static_cast<int>(ele_it - atom_data::element.begin()) + 1;
                        (*atom_it).m = 1; // Atom_Data::atomic_mass[(*atom_it).atomic_number - 1];
                        break;
                    }
                }
            }
            mol_ctr++;
        }
    }

    sys.num_part = 0;
    for (auto mi = 0; mi < sys.molecules.size(); mi++)
        sys.num_part += sys.molecules[mi].atoms.size();

}

void process::print() {
    std::cout << "\nNormal termination. Congratulations!" << std::endl;
}

void process::close() {
    // close files
    in.close();
    out.close();
}