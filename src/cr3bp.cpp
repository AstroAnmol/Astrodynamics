#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include "cr3bp.h"

//cr3bp::cr3bp(){}

void cr3bp::set_init_rot_coor(Eigen::Vector3d r, Eigen::Vector3d v){
    R=r;
    V=v;
}

void cr3bp::print_init_rot_coor(){
    std::cout<< "Radius:"<<R.transpose()<< std::endl;
    std::cout<< "Velocity:"<<V.transpose()<< std::endl;
}

void cr3bp::set_mu(double a){
    mu=a;
}

void cr3bp::set_DU(double a){
    DU=a;
}

Eigen::VectorXd cr3bp::EoM_rot(Eigen::Vector3d r, Eigen::Vector3d v){
    Eigen::VectorXd acc(6);
    Eigen::Vector3d R1, R2, R13, R23;
    R1<<-mu,0,0;
    R2<<1-mu,0,0;
    R13=-R1+r;
    R23=-R2+r;
    acc(0)=v(0);
    acc(1)=v(1);
    acc(2)=v(2);
    acc(3)=r(0) + 2*v(1) - (1-mu)*(r(0)-R1(0))/std::pow(R13.norm(),3) - mu*(r(0)-R2(0))/std::pow(R23.norm(),3);
    acc(4)=-2*v(0) + r(1) - (1-mu)*r(1)/std::pow(R13.norm(),3) - mu*r(1)/std::pow(R23.norm(),3);
    acc(5)=-(1-mu)*r(2)/std::pow(R13.norm(),3) - mu*r(2)/std::pow(R23.norm(),3);
    return acc;
}

void cr3bp::propagate(double step, double trange, std::string name){
    int t=0;
    //number of steps
    int noi=trange/step;

    Eigen::ArrayXXd R_propagated(3,noi+1);
    Eigen::ArrayXXd V_propagated(3,noi+1);
    Eigen::ArrayXXd R_inertial(3,noi+1);
    Eigen::ArrayXXd R_primary(3,noi+1);
    Eigen::ArrayXXd R_secondary(3,noi+1);

    Eigen::ArrayXXd time(1,noi+1);

    Eigen::ArrayXXd y(6,1);
    Eigen::ArrayXXd k1(6,1);
    Eigen::ArrayXXd k2(6,1);
    Eigen::ArrayXXd k3(6,1);
    Eigen::ArrayXXd k4(6,1);

    Eigen::Matrix3d TRI;

    Eigen::Vector3d RP, RS;
    RP<<-mu,0,0;
    RS<<1-mu,0,0;

    // setting initial values
    R_propagated.col(0)=R.transpose();
    V_propagated.col(0)=V.transpose();
    R_inertial.col(0)=R.transpose();
    R_primary.col(0)=RP.transpose();
    R_secondary.col(0)=RS.transpose();
    time.col(0)=0;
    y.block(0,0,3,1)=R;
    y.block(3,0,3,1)=V;
    while (t<noi){
        //RK4 Integrator
        k1= step*EoM_rot(R_propagated.col(t),V_propagated.col(t));
        k2= step*EoM_rot(R_propagated.col(t) + k1.block(0,0,3,1)/2,V_propagated.col(t)+ k1.block(3,0,3,1)/2);
        k3= step*EoM_rot(R_propagated.col(t) + k2.block(0,0,3,1)/2,V_propagated.col(t)+ k2.block(3,0,3,1)/2);
        k4= step*EoM_rot(R_propagated.col(t) + k3.block(0,0,3,1),V_propagated.col(t)+ k3.block(3,0,3,1));
        y += (k1+2*k2+2*k3+k4)/6;

        //Increasing Time Step
        t++;
        //set time
        time.col(t)=step*t;

        TRI<<   std::cos(step*t),std::sin(step*t),0,
                -std::sin(step*t),std::cos(step*t),0,
                0,0,1;
        //set R, V and Acc
        R_propagated.col(t)=y.block(0,0,3,1);
        V_propagated.col(t)=y.block(3,0,3,1);
        Eigen::Vector3d R_rot;
        R_rot=R_propagated.col(t);
        R_inertial.col(t)=TRI*R_rot;
        R_primary.col(t)=TRI*RP;
        R_secondary.col(t)=TRI*RS;
    }
    
    // converting to km
    R_inertial=DU*R_inertial;
    R_primary=DU*R_primary;
    R_secondary=DU*R_secondary;

    // Eigen Matrix print format to write csv files
    Eigen::IOFormat csv(10, 0, ", ", "\n", "", "", "", "");
    // Write file with R vector;
    // Matrix of all data
    Eigen::ArrayXXd Matrice(noi+1,16);
    //time
    Matrice.col(0)=time.transpose();
    
    //Radius
    Matrice.col(1)=R_propagated.row(0);
    Matrice.col(2)=R_propagated.row(1);
    Matrice.col(3)=R_propagated.row(2);
    //Velocity
    Matrice.col(4)=V_propagated.row(0);
    Matrice.col(5)=V_propagated.row(1);
    Matrice.col(6)=V_propagated.row(2);
    //Inertial Radius
    Matrice.col(7)=R_inertial.row(0);
    Matrice.col(8)=R_inertial.row(1);
    Matrice.col(9)=R_inertial.row(2);
    //Primary Body Radius
    Matrice.col(10)=R_primary.row(0);
    Matrice.col(11)=R_primary.row(1);
    Matrice.col(12)=R_primary.row(2);
    //Secondary Body Radius
    Matrice.col(13)=R_secondary.row(0);
    Matrice.col(14)=R_secondary.row(1);
    Matrice.col(15)=R_secondary.row(2);
    std::ofstream theFile;
    theFile.open(name + "_file.csv");
    theFile << "Time (TU),Radius_1 (DU),Radius_2 (DU),Radius_3 (DU),Velocity_1 (DU/TU),Velocity_2 (DU/TU),Velocity_3 (DU/TU),R_inertial_1 (km),R_inertial_2 (km),R_inertial_3 (km),R_prim_1 (km),R_prim_2 (km),R_prim_3 (km),R_sec_1 (km),R_sec_2 (km),R_sec_3 (km)" << std::endl;
    theFile << Matrice.format(csv)<< std::endl;
    theFile.close();
}