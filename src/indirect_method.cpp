#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include "indirect_method.h"
#include "orbit.h"

//default constructor
/* IndirectMethod::IndirectMethod(){

} */
//constructor
IndirectMethod::IndirectMethod(Eigen::VectorXd param){
    //GAparam=param;
    t0=param(0);
    m0=param(1);
    double Vinf_0_rad, Vinf_0_east, Vinf_0_north;
    Vinf_0_rad=param(2);
    Vinf_0_north=param(3);
    Lambda_R0(0)=param(4);
    Lambda_R0(1)=param(5);
    Lambda_R0(2)=param(6);
    tf=param(7);

    double muE=398600;//km^3/s^2
    double muS=132712440018;//km^3/s^2


    double AU=149597870.700; //km
    double TU=std::sqrt(AU*AU*AU/(muE+muS));

    //Setting Mars orbit
    // Initial orbit parameters taken from JPL HORIZONS in KM/S
    RM0<< 2.069270543147017e+08, -3.560689745239088E+06, -5.147936537447235e+06;
    VM0<< 1.304308833322233e+00, 2.628158890420931e+01, 5.188465740839767e-01;
    Mars.set_mu(1);
    Mars.set_cartesian(RM0, VM0);
    Eigen::VectorXd RVM(6);
    RVM=Mars.Kepler_prob(RM0,VM0,tf*TU);
    RMtf=RVM.block(0,0,3,1);
    VMtf=RVM.block(3,0,3,1);

    //Setting Earth orbit
    // Initial orbit parameters taken from JPL HORIZONS in KM/S
    RE0<<-2.627892928682480e+07, 1.445102393586391e+08, 3.022818135935813e+04;
    VE0<<-2.983052803283506e+01, -5.220465685407924e+00, -1.014621798034465e-04;
    Earth.set_mu(1);
    Earth.set_cartesian(RE0,VE0);
    Eigen::VectorXd RVE(6);
    RVE=Earth.Kepler_prob(RE0,VE0,t0*TU);
    REt0=RVE.block(0,0,3,1);
    VEt0=RVE.block(3,0,3,1);

    vc=std::sqrt(398600/(6371+200));// km/s
    ve=std::sqrt(2*398600/(6371+200));// km/s
    
    // non-dimentionalising the vectors
    vc=vc*TU/AU;
    ve=ve*TU/AU;
    RMtf=RMtf/AU;
    VMtf=VMtf*TU/AU;
    REt0=REt0/AU;
    VEt0=VEt0*TU/AU;

    mLEO=1/((1+epsilon)*std::exp((vc-ve)/c_dash)-epsilon);
    del_v=-c_dash*std::log((1/(1+epsilon))*(m0/mLEO + epsilon));

    v_inf_0=std::sqrt(std::pow(del_v + vc, 2) - std::pow(ve,2));
    Vinf_0_east=std::sqrt(v_inf_0*v_inf_0 - Vinf_0_north*Vinf_0_north - Vinf_0_rad*Vinf_0_rad);
    V_inf_0=Vinf_0_rad*REt0/REt0.norm() + Vinf_0_east*VEt0/VEt0.norm() + Vinf_0_north*(VEt0.cross(REt0))/(VEt0.norm()*REt0.norm());

    R0=REt0;
    V0=V_inf_0 + VEt0;

    Lambda_V0=V_inf_0*(m0+epsilon)*lambda_m0/(c_dash*(std::sqrt(std::pow(V_inf_0.norm(),2)+std::pow(ve,2))));

    Eigen::VectorXd x0;
    x0=Eigen::VectorXd::Zero(14);


    x0.block(0,0,3,1)=R0;
    x0.block(3,0,3,1)=V0;
    x0(6)=m0;
    x0.block(7,0,3,1)=Lambda_R0;
    x0.block(10,0,3,1)=Lambda_V0;
    x0(13)=lambda_m0;

    propagate();
}

double IndirectMethod::getFitness(){
    return fitness;
}

Eigen::Vector3d IndirectMethod::Thrust(Eigen::VectorXd x){//x is R,V,m,Lambda_R, Lambda_V, lambda_m
    Eigen::Vector3d R, V, Lambda_R, Lambda_V;
    double m, lambda_m;

    R=x.block(0,0,3,1);
    V=x.block(3,0,3,1);
    m=x(6);
    Lambda_R=x.block(7,0,3,1);
    Lambda_V=x.block(10,0,3,1);
    lambda_m=x(13);

    double r, lambda_V, Tmax;
    r=R.norm();
    lambda_V=Lambda_V.norm();

    Tmax=T0/(r*r);

    double SF;
    SF=lambda_V/m -lambda_m/c;

    Eigen::Vector3d T_out;

    if (SF>0){
        T_out=Tmax*Lambda_V/lambda_V;
    }
    else if (SF<0){
        T_out<<0,0,0;
    }
    return T_out;
}

Eigen::VectorXd IndirectMethod::differential(Eigen::VectorXd x){//x is R,V,m,Lambda_R, Lambda_V, lambda_m
    Eigen::Vector3d R, V, Lambda_R, Lambda_V;
    double m, lambda_m;

    R=x.block(0,0,3,1);
    V=x.block(3,0,3,1);
    m=x(6);
    Lambda_R=x.block(7,0,3,1);
    Lambda_V=x.block(10,0,3,1);
    lambda_m=x(13);

    double r, lambda_V;
    r=R.norm();
    lambda_V=Lambda_V.norm();

    Eigen::Vector3d T;
    T=Thrust(x);
    double t;
    t=T.norm();

    Eigen::VectorXd output(14);

    output(0)=V(0);
    output(1)=V(1);
    output(2)=V(2);
    output(3)=-R(0)/std::pow(r,3) + T(0)/m;
    output(4)=-R(1)/std::pow(r,3) + T(1)/m;
    output(5)=-R(2)/std::pow(r,3) + T(2)/m;
    output(6)=-t/c;
    output(7)=-Lambda_V(0)/std::pow(r,3) + 3*R(0)*(R.dot(Lambda_V))/std::pow(r,5);
    output(8)=-Lambda_V(1)/std::pow(r,3) + 3*R(1)*(R.dot(Lambda_V))/std::pow(r,5);
    output(9)=-Lambda_V(2)/std::pow(r,3) + 3*R(2)*(R.dot(Lambda_V))/std::pow(r,5);
    output(10)=-Lambda_R(0);
    output(11)=-Lambda_R(1);
    output(12)=-Lambda_R(2);
    output(13)=t*lambda_V/(m*m);
    return output;
}

void IndirectMethod::propagate(){
    int t=0;
    double step=0.01;
    //number of steps
    int noi=tf/step;

    Eigen::ArrayXXd x_propagated(14,noi+1);
    
    Eigen::ArrayXXd time(1,noi+1);

    Eigen::ArrayXXd k1(14,1);
    Eigen::ArrayXXd k2(14,1);
    Eigen::ArrayXXd k3(14,1);
    Eigen::ArrayXXd k4(14,1);

    Eigen::VectorXd x0;
    x0=Eigen::VectorXd::Zero(14);

    x0.block(0,0,3,1)=R0;
    x0.block(3,0,3,1)=V0;
    x0(6)=m0;
    x0.block(7,0,3,1)=Lambda_R0;
    x0.block(10,0,3,1)=Lambda_V0;
    x0(13)=lambda_m0;

    // setting initial values
    x_propagated.col(0)=x0;
    time.col(0)=0;
    while (t<noi){
        //RK4 Integrator
        k1= step*differential(x_propagated.col(t));
        k2= step*differential(x_propagated.col(t) + k1/2);
        k3= step*differential(x_propagated.col(t) + k2/2);
        k4= step*differential(x_propagated.col(t) + k3);
        x_propagated.col(t+1)= x_propagated.col(t)+(k1+2*k2+2*k3+k4)/6;

        //Increasing Time Step
        t++;
        //set time
        time.col(t)=step*t;
    }
    Eigen::VectorXd xf;
    Eigen::Vector3d delR_f, delV_f;
    xf=x_propagated.col(t);
    delR_f=(xf.block(0,0,3,1)-RMtf);
    delV_f=(xf.block(3,0,3,1)-VMtf);

    delR_f=delR_f.cwiseAbs();
    delV_f=delV_f.cwiseAbs();

    fitness=delR_f.sum()+delV_f.sum();
}

void IndirectMethod::save(std::string name){
    int t=0;
    double step=0.01;
    //number of steps
    int noi=tf/step;

    Eigen::ArrayXXd x_propagated(14,noi+1);
    
    Eigen::ArrayXXd time(1,noi+1);

    Eigen::ArrayXXd k1(14,1);
    Eigen::ArrayXXd k2(14,1);
    Eigen::ArrayXXd k3(14,1);
    Eigen::ArrayXXd k4(14,1);

    Eigen::VectorXd x0;
    x0=Eigen::VectorXd::Zero(14);

    x0.block(0,0,3,1)=R0;
    x0.block(3,0,3,1)=V0;
    x0(6)=m0;
    x0.block(7,0,3,1)=Lambda_R0;
    x0.block(10,0,3,1)=Lambda_V0;
    x0(13)=lambda_m0;

    // setting initial values
    x_propagated.col(0)=x0;
    time.col(0)=0;
    while (t<noi){
        //RK4 Integrator
        k1= step*differential(x_propagated.col(t));
        k2= step*differential(x_propagated.col(t) + k1/2);
        k3= step*differential(x_propagated.col(t) + k2/2);
        k4= step*differential(x_propagated.col(t) + k3);
        x_propagated.col(t+1)= x_propagated.col(t)+(k1+2*k2+2*k3+k4)/6;

        //Increasing Time Step
        t++;
        //set time
        time.col(t)=step*t;
    }
    // Eigen Matrix print format to write csv files
    Eigen::IOFormat csv(10, 0, ", ", "\n", "", "", "", "");
    // Write file with R vector;
    // Matrix of all data
    Eigen::ArrayXXd Matrice(noi+1,14);
    Matrice=x_propagated.transpose();
    Matrice.conservativeResize(Matrice.rows(),Matrice.cols()+1);
    Matrice.col(Matrice.cols()-1)=time.transpose();
    std::ofstream theFile;
    theFile.open(name + "_file.csv");
    theFile << "Earth Radius:," << REt0.format(csv) << std::endl;
    theFile << "Earth Velocity:," << VEt0.format(csv) << std::endl;
    theFile << "Radius_1, Radius_2, Radius_3, Velocity_1, Velocity_2, Velocity_3, Mass, Adjoint Variables..." << std::endl;
    theFile << Matrice.format(csv)<< std::endl<< std::endl<< std::endl;
    theFile << "Mars Radius:," << RMtf.format(csv) << std::endl;
    theFile << "Mars Velocity:," << VMtf.format(csv) << std::endl;
    theFile.close();
}