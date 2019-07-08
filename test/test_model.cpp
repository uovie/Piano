#include <iostream>
#include <fstream>

#include <Eigen/Dense>

#include "model.h"

int main()
{
    std::ofstream out0, out1, out2;
    out0.open("r.dat");
    out1.open("LJV.dat");
    out2.open("LJF.dat");

    Eigen::ArrayXXd r1 = Eigen::ArrayXXd::Zero(3, 1);
    Eigen::ArrayXXd r2 = Eigen::ArrayXXd::Zero(3, 1);

    r2(0, 0) = 2.6;

    uovie::model::Lennard_Jones LJ(35.6, 2.749);
    for (int i = 0; i < 100; i++) {
        out0 << r2(0, 0) << std::endl;
        out1 << LJ.V(r1, r2) << std::endl;
        out2 << (LJ.F(r1, r2))(0, 0) << std::endl;

        r2(0, 0) += 0.05;
    }

    return 0;
}