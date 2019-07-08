#include "model.h"

namespace uovie {
namespace model {

    double Lennard_Jones_Gaussian::V(const Eigen::ArrayXXd& ri, const Eigen::ArrayXXd& rj) {
        double pot = 0;
        for (int i = 0; i < c.size(); i++)
            pot += c[i] * exp(-alpha[i] * ((ri - rj).pow(2)).sum() / pow(sigma, 2));
        pot *= epsilon;

        return pot;
    }

    Eigen::ArrayXXd Lennard_Jones_Gaussian::F(const Eigen::ArrayXXd& ri, const Eigen::ArrayXXd& rj) {
        double tmp_term = 0;
        for (int i = 0; i < c.size(); i++)
            tmp_term -= 2 * c[i] * alpha[i] * pow(sigma, -2) * exp(-alpha[i] * ((ri - rj).pow(2)).sum() / pow(sigma, 2));
        tmp_term *= epsilon;

        return tmp_term * (ri - rj);
    }
    
} // !model
} // !uovie