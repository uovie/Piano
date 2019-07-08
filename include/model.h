/* Harmonic Bead-Spring Model */
#ifndef MODEL_H_
#define MODEL_H_

#include <vector>
#include <cmath>
#include <Eigen/Dense>

namespace uovie {
namespace model {

    class harmonic_oscilator {
    public:
        harmonic_oscilator() = default;
        harmonic_oscilator(const double& _omega): omega(_omega) { }

        double ome() const { return omega; }

    private:
        const double omega;
    };

    class Lennard_Jones {
    public:
        Lennard_Jones() = default;
        Lennard_Jones(const double& _epsilon, const double& _sigma):
            epsilon(_epsilon), sigma(_sigma) { }

        double eps() const { return epsilon; }
        double sig() const { return sigma; }
        double V(const Eigen::ArrayXXd& ri, const Eigen::ArrayXXd& rj) {
            return 4 * epsilon * (pow(sigma, 12) / pow(((ri - rj).pow(2)).sum(), 6)
                - pow(sigma, 6) / pow(((ri - rj).pow(2)).sum(), 3));
        }
        Eigen::ArrayXXd F(const Eigen::ArrayXXd& ri, const Eigen::ArrayXXd& rj) {
            return 24 * epsilon * (2 * pow(sigma, 12) / pow(((ri - rj).pow(2)).sum(), 7)
                - pow(sigma, 6) / pow(((ri - rj).pow(2)).sum(), 4)) * (ri - rj);
        }

    private:
        const double epsilon;
        const double sigma;
    };

    class Lennard_Jones_Gaussian {
    public:
        Lennard_Jones_Gaussian() = default;
        Lennard_Jones_Gaussian(const double& _epsilon, const double& _sigma,
            const std::vector<double>& _c, const std::vector<double>& _alpha) :
            epsilon(_epsilon), sigma(_sigma), c(_c), alpha(_alpha) { }

        double eps() const { return epsilon; }
        double sig() const { return sigma; }
        double coe(int& i) const { return c[i]; }
        double alp(int& i) const { return alpha[i]; }
        double V(const Eigen::ArrayXXd& ri, const Eigen::ArrayXXd& rj);
        Eigen::ArrayXXd F(const Eigen::ArrayXXd& ri, const Eigen::ArrayXXd& rj);
         
    private:
        const double epsilon;
        const double sigma;
        const std::vector<double> c;
        const std::vector<double> alpha;
    };

} // !model
} // !uovie

#endif // !HARMONIC_BEAD_SPRING_H_