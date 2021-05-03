#ifndef Genetic_GA_H
#define Genetic_GA_H
#include <eigen-3.3.7/Eigen/Dense>

class Genetic_GA{
    public:
        Genetic_GA(std::string seq);
        double get_BestFitness();
        Eigen::VectorXd get_BestGene();
    private:
        int noi=50; // number of iterations
        int pop_size=100; // population size
        double t0_lb=58;
        double t0_up=59;
        double tf_lb=73;
        double tf_up=74;
        double m0_lb=0.97;
        double m0_up=1;
        double lambdaR0_lb=-2;
        double lambdaR0_ub=2;
        double vinf_lp=-0.025;//8.3936E-4;
        double vinf_up=0.025;//8.3936E-4;

        double BestFitness;
        Eigen::VectorXd BestGene;

        Eigen::Vector3i top3(Eigen::VectorXd a);
        double fRand(double fMin, double fMax);
        double beta_dash(double u);
};
#endif