#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <eigen-5.0.0/Eigen/Dense>
#include <cmath>
#include <queue>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include "orbit.h"
#include "debris.h"
// #include "cr3bp.h"
// #include "Indirect_BVP_DM.h"
// #include "Genetic_DM.h"
// #include "Indirect_BVP_GA.h"
// #include "Genetic_GA.h"

int main(){
        // srand(time(0));
        /*
        Eigen::VectorXd Gene(8);
        Gene<<50.6555, 0.972694, -0.000650264, -0.000551347, -0.590737, -1.39754, -0.552314, 63.0691;
        IndirectMethod Trial(Gene);
        Trial.print("VEt0");
        Trial.print("Vinf");
        Trial.print("all scalars");
        Trial.save("Try");
        */
        /*
        for (int i = 0; i < 1; i++)
        {
            Genetic_GA trial1("EEMM");
            Eigen::VectorXd BestGeneT1(24);
            //BestGeneT1<<66, 0.972694, -0.000650264, -0.000551347, -0.590737, -1.39754, -0.552314, 74, -0.000650264, -0.000551347, -0.590737, -1.39754, -0.552314, 79, -0.000650264, -0.000551347, -0.590737, -1.39754, -0.552314, -0.590737, -1.39754, -0.552314, 95, 72;
            BestGeneT1=trial1.get_BestGene();
            std::cout<<BestGeneT1.transpose()<<std::endl;
            std::cout<<trial1.get_BestFitness()<<std::endl;
            Indirect_BVP_GA Best(BestGeneT1, "EEMM");
            std::string name;
            std::ostringstream oss;
            oss << "EEMM_" << i;
            name=oss.str();
            Best.save(name);
        }    
        */
    // satellite orbit
    double a, e, i, omega, Omega, theta;
    a =         750 + 6371; // km;
    e =         0.063;
    i =         45;
    omega =    0;
    Omega =     0;
    theta =     0;

    orbit sat;
    Eigen::VectorXd OE(6);
    OE << a, e, i, omega, Omega, theta;
    sat.set_OE(OE);
    sat.set_mu(0);
    std::cout<< "Satellite initial Cartesian state: \n";
    sat.print_cartesian();
    double time_period = sat.get_TimePeriod();
    // propagate satellite orbit at higher time step for plot
    double step_sat = 10; // seconds
    // sat.propagate_2BP(step_sat, time_period, 0, "sat_orbit");

    // soliton characteristics
    double cone_angle = 45 * M_PI / 180; // radians
    double cone_height = 10;              // km

    // debris object
    Eigen::Vector3d r1, v1;
    double a_d, e_d, i_d, omega_d, Omega_d, theta_d;
    a_d =         a; // km;
    e_d =         e;
    i_d =         180-i;
    omega_d =     180 + omega;
    Omega_d =     180 + Omega;
    theta_d =     theta - 0.01;

    Eigen::VectorXd OE_d(6);
    OE_d << a_d, e_d, i_d, omega_d, Omega_d, theta_d;
    orbit debris;
    debris.set_OE(OE_d);
    debris.set_mu(0);
    debris.print_cartesian();
    Eigen::VectorXd cartesian_d = debris.get_cartesian();
    r1 = cartesian_d.segment(0,3);
    v1 = cartesian_d.segment(3,3);
    Debris D1;
    D1.set_state(r1, v1);
    D1.set_soliton_cone(cone_angle, cone_height);

    std::cout << "Soliton velocity: " << D1.get_soliton_vel() << " km/s" << std::endl;
    // get time to reach cone base
    double time_to_cone_base = D1.get_time_to_reach_cone_base();
    std::cout << "Time to reach cone base: " << time_to_cone_base << " seconds" << std::endl;


    // propagate satellite orbit to the same time
    double step_size = 0.001; // seconds
    sat.propagate_2BP(step_size, time_to_cone_base, 0, "sat_prop");


    D1.read_sat_orbit("sat_prop");

    D1.check_detection();

}