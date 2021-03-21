#ifndef CR3BP_H
#define CR3BP_H
#include <eigen-3.3.7/Eigen/Dense>

class cr3bp{
    private:
        //Rotating Frame Coordinates
		Eigen::Vector3d R;
		Eigen::Vector3d V;
        double mu, DU;
        Eigen::VectorXd EoM_rot(Eigen::Vector3d r, Eigen::Vector3d v);
    public:
        //cr3bp(/* args */);
        void set_init_rot_coor(Eigen::Vector3d r, Eigen::Vector3d v);
        void set_mu(double a);
        void set_DU(double a);
        //void get_init_rot_coor();
        void print_init_rot_coor();

        void propagate(double step, double trange, std::string name);
};

#endif