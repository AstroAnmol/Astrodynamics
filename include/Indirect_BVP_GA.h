#ifndef IND_GA_H
#define IND_GA_H
#include <eigen-3.3.7/Eigen/Dense>
#include "orbit.h"

class Indirect_BVP_GA{
    public:
        /* Indirect_BVP_GA(); */
        Indirect_BVP_GA(Eigen::VectorXd param, std::string seq);
        double getFitness();
        Eigen::Vector3d Thrust(Eigen::VectorXd x);//x is R,V,m,Lambda_R, Lambda_V, lambda_m        
        Eigen::VectorXd differential(Eigen::VectorXd x);//x is R,V,m,Lambda_R, Lambda_V, lambda_m
        Eigen::VectorXd propagate(Eigen::VectorXd x0, double ti, double tf);//
        void propagate_save(Eigen::VectorXd x0, double ti, double tf, std::string name);//
        // to save the propagated trajectory
        void save(std::string name);
        // print some values
        void print(std::string value);
    private:
        // time instances of all four planetary encounters
        double t0, tGA1, tGA2, tf, tstar;
        // Spacecraft properties
        // Mass
        double m0, lambda_m0=1;
        double mGA1_i, lambda_mGA1_i;
        double mGA1_f, lambda_mGA1_f;
        double mGA2_i, lambda_mGA2_i;
        double mGA2_f, lambda_mGA2_f;
        double mf, lambda_mf;
        double mLEO;
        // State Vector and Adjoint Vector
        Eigen::Vector3d R0, V0, Lambda_R0, Lambda_V0;
        Eigen::Vector3d RGA1_i, VGA1_i, Lambda_RGA1_i, Lambda_VGA1_i;
        Eigen::Vector3d RGA1_f, VGA1_f, Lambda_RGA1_f, Lambda_VGA1_f;
        Eigen::Vector3d RGA2_i, VGA2_i, Lambda_RGA2_i, Lambda_VGA2_i;
        Eigen::Vector3d RGA2_f, VGA2_f, Lambda_RGA2_f, Lambda_VGA2_f;
        Eigen::Vector3d Rf, Vf, Lambda_Rf, Lambda_Vf;
        
        // stacked state and adjoint vectors
        Eigen::VectorXd x0, xGA1_i;
        Eigen::VectorXd xGA1_f, xGA2_i;
        Eigen::VectorXd xGA2_f, xf;

        // Escape Velocity (Origin 200 km orbit)
        // Circular orbit Velocity (Origin 200 km orbit)
        // Impusive Delta V to leave Earth orbit
        double ve, vc, del_v;

        // Circular velocity at GA2 (500 km orbit)
        double vp;

        // Hyperbolic Excess velocities        
        double v_inf_0;
        Eigen::Vector3d V_inf_0;
        double v_inf_GA1_i, v_inf_GA1_f;
        Eigen::Vector3d V_inf_GA1_i, V_inf_GA1_f;
        double v_inf_GA2_i, v_inf_GA2_f;
        Eigen::Vector3d V_inf_GA2_i, V_inf_GA2_f;

        // Propulsion properties
        const double T0=0.02; // Thrust at 1 AU
        const double c=1;   // exhaust velocity
        const double c_dash=0.0955;
        const double epsilon=0.06;

        // orbits of planets
        orbit Earth, Mars;

        //orbit destination
        Eigen::Vector3d RD0, VD0; // Destination vector on 01/01/2000
        orbit Destination;
        Eigen::Vector3d RDtf, VDtf;

        //Gravity Assist 1 destination
        Eigen::Vector3d RGA10, VGA10; // Destination vector on 01/01/2000
        orbit GA1;
        Eigen::Vector3d RGA1, VGA1;

        //Gravity Assist 2 destination
        Eigen::Vector3d RGA20, VGA20; // Destination vector on 01/01/2000
        orbit GA2;
        Eigen::Vector3d RGA2, VGA2;

        //orbit orgin
        Eigen::Vector3d RO0, VO0; // Origin vector on 01/01/2000
        orbit Origin;
        Eigen::Vector3d ROt0, VOt0;

        double fitness;
};


#endif