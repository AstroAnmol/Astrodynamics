#ifndef Genetic_DM_H
#define Genetic_DM_H
#include <eigen-3.3.7/Eigen/Dense>

class Genetic_DM{
    public:
        Genetic_DM();
        double get_BestFitness();
        Eigen::VectorXd get_BestGene();
    private:
        const int noi=50; // number of iterations
        const int pop_size=100; // population size
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
        
};


#endif