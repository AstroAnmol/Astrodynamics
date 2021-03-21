#ifndef GA_H
#define GA_H
#include <eigen-3.3.7/Eigen/Dense>
#include "orbit.h"
#include "indirect_method.h"

class GA{
    public:
        GA();
        double get_BestFitness();
        Eigen::VectorXd get_BestGene();
    private:
        const int noi=50; // number of iterations
        const int pop_size=100; // population size
        double t0_lb=50;
        double t0_up=63;
        double tf_lb=60;
        double tf_up=93;
        double m0_lb=0.97;
        double m0_up=1;
        double lambdaR0_lb=-2;
        double lambdaR0_ub=2;
        double vinf_lp=-0.025;
        double vinf_up=0.025;

        double BestFitness;
        Eigen::VectorXd BestGene;
        
};


#endif