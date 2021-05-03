#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include "Indirect_BVP_GA.h"
#include "orbit.h"

//default constructor
/* Indirect_BVP_GA::Indirect_BVP_GA(){

} */

// constructor takes GA paramter and Gravity assist seqeunce. For example: EEMM means Earth Dep, followed by Earth GA followed by Mars GA and finally Mars arrival.
Indirect_BVP_GA::Indirect_BVP_GA(Eigen::VectorXd param, std::string seq){
    // ------------------------------------------------------------
    // HAVE TO ADD A IF CONDITION FOR NUMBER OF PARAMETERS AND SEQ
    // ------------------------------------------------------------
    // set fitness to zero
    fitness=0;
    // Initial Origin Departure V Inf Vector
    double Vinf_0_rad, Vinf_0_east, Vinf_0_north;
    // First Gravity Assist Departure V Inf Vector
    double Vinf_GA1_f_rad, Vinf_GA1_f_east, Vinf_GA1_f_north;
    // Second Gravity Assist Departure V Inf Vector
    double Vinf_GA2_f_rad, Vinf_GA2_f_east, Vinf_GA2_f_north;
    //GAparam=param;
    {   // Origin/ Initial Parameters
        t0=param(0);
        m0=param(1);
        Vinf_0_rad=param(2);
        Vinf_0_north=param(3);
        Lambda_R0(0)=param(4);
        Lambda_R0(1)=param(5);
        Lambda_R0(2)=param(6);
        // First Gravity Assist Parameters
        tGA1=param(7);
        Vinf_GA1_f_rad=param(8);
        Vinf_GA1_f_north=param(9);
        Lambda_RGA1_f(0)=param(10);
        Lambda_RGA1_f(1)=param(11);
        Lambda_RGA1_f(2)=param(12);
        // Second Gravity Assist Parameters (with minimum height constraint)
        tGA2=param(13);
        Vinf_GA2_f_rad=param(14);
        Vinf_GA2_f_north=param(15);
        Lambda_RGA2_f(0)=param(16);
        Lambda_RGA2_f(1)=param(17);
        Lambda_RGA2_f(2)=param(18);
        Lambda_VGA2_f(0)=param(19);
        Lambda_VGA2_f(1)=param(20);
        Lambda_VGA2_f(2)=param(21);
        // Destination/ Final Parameters
        tf=param(22);
    }
    // Orbit Parameters of planets as on 2000-01-01
    {
        //Setting Earth orbit
        // Initial orbit parameters taken from JPL HORIZONS in KM/S
        Eigen::Vector3d RE0, VE0;
        RE0<<-2.521092855899356E+07, 1.449279195838006E+08, -6.164165719002485E+02;
        VE0<<-2.983983333677879E+01, -5.207633902410673E+00, 6.168441184239981E-05;
        Earth.set_mu(1);
        Earth.set_cartesian(RE0,VE0);

        //Setting Mars orbit
        // Initial orbit parameters taken from JPL HORIZONS in KM/S
        Eigen::Vector3d RM0, VM0;
        RM0<< 2.079950549908331E+08, -3.143009561106971E+06, -5.178781160069674E+06;
        VM0<< 1.295003532851602E+00, 2.629442067068712E+01, 5.190097267545717E-01;
        Mars.set_mu(1);
        Mars.set_cartesian(RM0, VM0);
    }
    //constants
    double muE=398600;//km^3/s^2
    double muM=42828; //km^3/s^2
    double muS=132712440018;//km^3/s^2
    double Re=6371; //km
    double Rm=3389.5; //km

    double AU=149597870.700; //km
    double VU=std::sqrt(muS/AU);
    double TU=AU/VU;
    
    //Setting Origin orbit (Earth for now, Later add based on Seq)
    // Initial orbit parameters taken from JPL HORIZONS in KM/S
    Origin=Earth;
    Eigen::VectorXd RVo0(6);
    RVo0=Origin.get_cartesian();
    RO0=RVo0.block(0,0,3,1);
    VO0=RVo0.block(3,0,3,1);
    Eigen::VectorXd RVo(6);
    RVo=Origin.Kepler_prob(RO0,VO0,t0*TU);
    ROt0=RVo.block(0,0,3,1);
    VOt0=RVo.block(3,0,3,1);
    double muO=muE;
    double Ro=Re;


    //Setting GA1 orbit (Earth for now, Later add based on Seq)
    // Initial orbit parameters taken from JPL HORIZONS in KM/S
    GA1=Earth;
    Eigen::VectorXd RVGA10(6);
    RVGA10=GA1.get_cartesian();
    RGA10=RVGA10.block(0,0,3,1);
    VGA10=RVGA10.block(3,0,3,1);
    Eigen::VectorXd RVGA1(6);
    RVGA1=GA1.Kepler_prob(RGA10,VGA10,tGA1*TU);
    RGA1=RVGA1.block(0,0,3,1);
    VGA1=RVGA1.block(3,0,3,1);

    //Setting GA2 orbit (Mars for now, Later add based on Seq)
    // Initial orbit parameters taken from JPL HORIZONS in KM/S
    GA2=Mars;
    Eigen::VectorXd RVGA20(6);
    RVGA20=GA2.get_cartesian();
    RGA20=RVGA20.block(0,0,3,1);
    VGA20=RVGA20.block(3,0,3,1);
    Eigen::VectorXd RVGA2(6);
    RVGA2=GA2.Kepler_prob(RGA20,VGA20,tGA2*TU);
    RGA2=RVGA2.block(0,0,3,1);
    VGA2=RVGA2.block(3,0,3,1);
    double muGA2=muM;
    double Rga2=Rm;

    //Setting Destination orbit (Mars for now, Later add based on Seq)
    // Initial orbit parameters taken from JPL HORIZONS in KM/S
    Destination=Mars;
    Eigen::VectorXd RVD0(6);
    RVD0=Destination.get_cartesian();
    RD0=RVD0.block(0,0,3,1);
    VD0=RVD0.block(3,0,3,1);
    Eigen::VectorXd RVD(6);
    RVD=Destination.Kepler_prob(RD0,VD0,tf*TU);
    RDtf=RVD.block(0,0,3,1);
    VDtf=RVD.block(3,0,3,1);

    ve=std::sqrt(2*muO/(Ro+200));// km/s
    vc=std::sqrt(muO/(Ro+200));// km/s
    
    vp=std::sqrt(2*muGA2/(Rga2+500));// km/s

    // non-dimentionalising the vectors
    vc =    vc/VU;
    ve =    ve/VU;
    vp =    vp/VU;
    RDtf =  RDtf/AU;
    VDtf =  VDtf/VU;
    ROt0 =  ROt0/AU;
    VOt0 =  VOt0/VU;
    RGA1 =  RGA1/AU;
    VGA1 =  VGA1/VU;

    // Initial Boundary Conditions
    mLEO=1/((1+epsilon)*std::exp((vc-ve)/c_dash)-epsilon);
    del_v=-c_dash*std::log((1/(1+epsilon))*(m0/mLEO + epsilon));

    v_inf_0=std::sqrt(std::pow(del_v + vc, 2) - std::pow(ve,2));
    Vinf_0_east=std::sqrt(v_inf_0*v_inf_0 - Vinf_0_north*Vinf_0_north - Vinf_0_rad*Vinf_0_rad);
    V_inf_0=Vinf_0_rad*ROt0/ROt0.norm() + Vinf_0_east*VOt0/VOt0.norm() + Vinf_0_north*(VOt0.cross(ROt0))/(VOt0.norm()*ROt0.norm());

    R0=ROt0;
    V0=V_inf_0 + VOt0;

    Lambda_V0=V_inf_0*(m0+epsilon)*lambda_m0/(c_dash*(std::sqrt(std::pow(V_inf_0.norm(),2)+std::pow(ve,2))));

    // stacked state and adjoint vector
    
    x0=Eigen::VectorXd::Zero(14);

    x0.block(0,0,3,1)=R0;
    x0.block(3,0,3,1)=V0;
    x0(6)=m0;
    x0.block(7,0,3,1)=Lambda_R0;
    x0.block(10,0,3,1)=Lambda_V0;
    x0(13)=lambda_m0;
    // propagate first section of heliocentric trajectory
    xGA1_i=propagate(x0, t0, tGA1);
    // Define spacecraft vectors just before gravity assist
    RGA1_i=xGA1_i.block(0,0,3,1);
    VGA1_i=xGA1_i.block(3,0,3,1);
    mGA1_i=xGA1_i(6);
    Lambda_RGA1_i=xGA1_i.block(7,0,3,1);
    Lambda_VGA1_i=xGA1_i.block(10,0,3,1);
    lambda_mGA1_i=xGA1_i(13);

    // Check for constraints
        // RGA1_i==RGA1
        Eigen::Vector3d delR_GA1_i;
        delR_GA1_i=RGA1_i-RGA1;
        delR_GA1_i=delR_GA1_i.cwiseAbs();
        fitness=fitness+delR_GA1_i.sum();
    
    // Apply Boundary conditions for next section of Heliocentric Trajectory
    // Defining Position vectors
        RGA1_f=RGA1;
        // adjoint variable already defined at the start.
    // Defining Velocity vectors
        //initial hyperbolic excess velocity
        V_inf_GA1_i=VGA1_i-VGA1; //vector
        v_inf_GA1_i=V_inf_GA1_i.norm(); //magnitude
        // defining final hyperbolic excess velocity
        v_inf_GA1_f=v_inf_GA1_i; //magnitude
        Vinf_GA1_f_east=std::sqrt(v_inf_GA1_f*v_inf_GA1_f - Vinf_GA1_f_north*Vinf_GA1_f_north - Vinf_GA1_f_rad*Vinf_GA1_f_rad);
        V_inf_GA1_f=Vinf_GA1_f_rad*RGA1/RGA1.norm() + Vinf_GA1_f_east*VGA1/VGA1.norm() + Vinf_GA1_f_north*(VGA1.cross(RGA1))/(VGA1.norm()*RGA1.norm());
        //final spacecraft heliocentric velocity
        VGA1_f=V_inf_GA1_f+VGA1;
        // adjoint variable
        Lambda_VGA1_f=Lambda_VGA1_i.norm()*V_inf_GA1_f/v_inf_GA1_f;
    // Defining Mass
        mGA1_f=mGA1_i;
        lambda_mGA1_f=lambda_mGA1_i;
    
    // stacked state and adjoint vector
    xGA1_f=Eigen::VectorXd::Zero(14);
    xGA1_f.block(0,0,3,1)=RGA1_f;
    xGA1_f.block(3,0,3,1)=VGA1_f;
    xGA1_f(6)=mGA1_f;
    xGA1_f.block(7,0,3,1)=Lambda_RGA1_f;
    xGA1_f.block(10,0,3,1)=Lambda_VGA1_f;
    xGA1_f(13)=lambda_mGA1_f;
    // propagate second section of heliocentric trajectory
    xGA2_i=propagate(xGA1_f, tGA1, tGA2);
    // Define spacecraft vectors just before gravity assist
    RGA2_i=xGA2_i.block(0,0,3,1);
    VGA2_i=xGA2_i.block(3,0,3,1);
    mGA2_i=xGA2_i(6);
    Lambda_RGA2_i=xGA2_i.block(7,0,3,1);
    Lambda_VGA2_i=xGA2_i.block(10,0,3,1);
    lambda_mGA2_i=xGA2_i(13);

    // Check for position contraints
        // RGA2_i==RGA2
        Eigen::Vector3d delR_GA2_i;
        delR_GA2_i=RGA2_i-RGA2;
        delR_GA2_i=delR_GA2_i.cwiseAbs();
        fitness=fitness+delR_GA2_i.sum();

    // Apply Boundary conditions for next section of Heliocentric Trajectory
    // Defining position vectors
        RGA2_f=RGA2;
        // adjoint variable already defined at the start.
    // Defining Velocity vectors
        //initial hyperbolic excess velocity
        V_inf_GA2_i=VGA2_i-VGA2; //vector
        v_inf_GA2_i=V_inf_GA2_i.norm(); //magnitude
        // defining final hyperbolic excess velocity
        v_inf_GA2_f=v_inf_GA2_i; //magnitude
        Vinf_GA2_f_east=std::sqrt(v_inf_GA2_f*v_inf_GA2_f - Vinf_GA2_f_north*Vinf_GA2_f_north - Vinf_GA2_f_rad*Vinf_GA2_f_rad);
        V_inf_GA2_f=Vinf_GA2_f_rad*RGA2/RGA2.norm() + Vinf_GA2_f_east*VGA2/VGA2.norm() + Vinf_GA2_f_north*(VGA2.cross(RGA2))/(VGA2.norm()*RGA2.norm());
        //final spacecraft heliocentric velocity
        VGA2_f=V_inf_GA2_f+VGA2;
        // adjoint variable defined earlier
    // Defining Mass
        mGA2_f=mGA2_i;
        lambda_mGA2_f=lambda_mGA2_i;

    // Check for turn angle constraint
        // Vinf_f.dot(Vinf_i)=-cos(2*phi)*vinf_i^2
        // cos(phi)=vp^2/(vinf_i^2 + vp^2)
        double cos_phi, cos_2phi;
        cos_phi=vp*vp/(v_inf_GA2_i*v_inf_GA2_i + vp*vp);
        cos_2phi=2*cos_phi*cos_phi - 1;
        double error;
        error= V_inf_GA2_f.dot(V_inf_GA2_i) + cos_2phi*v_inf_GA2_i*v_inf_GA2_i;
        fitness= fitness + abs(error);
    
    // stacked state and adjoint vector
    xGA2_f=Eigen::VectorXd::Zero(14);
    xGA2_f.block(0,0,3,1)=RGA2_f;
    xGA2_f.block(3,0,3,1)=VGA2_f;
    xGA2_f(6)=mGA2_f;
    xGA2_f.block(7,0,3,1)=Lambda_RGA2_f;
    xGA2_f.block(10,0,3,1)=Lambda_VGA2_f;
    xGA2_f(13)=lambda_mGA2_f;
    // propagate second section of heliocentric trajectory
    xf=propagate(xGA2_f, tGA2, tf);
    // Define spacecraft vectors just before gravity assist
    Rf=xf.block(0,0,3,1);
    Vf=xf.block(3,0,3,1);
    mf=xf(6);
    Lambda_Rf=xf.block(7,0,3,1);
    Lambda_Vf=xf.block(10,0,3,1);
    lambda_mf=xf(13);
    
    // Check final constraints
        Eigen::Vector3d delR_f, delV_f;
        delR_f=(xf.block(0,0,3,1)-RDtf);
        delV_f=(xf.block(3,0,3,1)-VDtf);

        delR_f=delR_f.cwiseAbs();
        delV_f=delV_f.cwiseAbs();

        fitness=fitness + delR_f.sum()+delV_f.sum();

    if (std::isnan(fitness)==1)
    {
        fitness=1000;
    }
}

double Indirect_BVP_GA::getFitness(){
    return fitness;
}

Eigen::Vector3d Indirect_BVP_GA::Thrust(Eigen::VectorXd x){//x is R,V,m,Lambda_R, Lambda_V, lambda_m
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

Eigen::VectorXd Indirect_BVP_GA::differential(Eigen::VectorXd x){//x is R,V,m,Lambda_R, Lambda_V, lambda_m
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
    output(7)=(Lambda_V(0)/std::pow(r,3) - 3*R(0)*(R.dot(Lambda_V))/std::pow(r,5));
    output(8)=(Lambda_V(1)/std::pow(r,3) - 3*R(1)*(R.dot(Lambda_V))/std::pow(r,5));
    output(9)=(Lambda_V(2)/std::pow(r,3) - 3*R(2)*(R.dot(Lambda_V))/std::pow(r,5));
    output(10)=-Lambda_R(0);
    output(11)=-Lambda_R(1);
    output(12)=-Lambda_R(2);
    output(13)=t*lambda_V/(m*m);
    return output;
}
// propagate initial vector x0 from initial time ti to final time tf.
Eigen::VectorXd Indirect_BVP_GA::propagate(Eigen::VectorXd x0, double ti, double tf){
    int t=0;
    double step=0.01;
    //number of steps
    int noi=(tf-ti)/step;

    Eigen::ArrayXXd x_propagated(14,noi+1);
    
    Eigen::ArrayXXd time(1,noi+1);

    Eigen::ArrayXXd k1(14,1);
    Eigen::ArrayXXd k2(14,1);
    Eigen::ArrayXXd k3(14,1);
    Eigen::ArrayXXd k4(14,1);

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
    xf=x_propagated.col(t);
    return xf;
}

// propagate initial vector x0 from initial time ti to final time tf for saving.
void Indirect_BVP_GA::propagate_save(Eigen::VectorXd x0, double ti, double tf, std::string name){
    int t=0;
    double step=0.01;
    //number of steps
    int noi=(tf-ti)/step;

    Eigen::ArrayXXd x_propagated(14,noi+1);
    
    
    Eigen::ArrayXXd time(1,noi+1);
    Eigen::ArrayXXd th(1,noi+1);

    Eigen::ArrayXXd k1(14,1);
    Eigen::ArrayXXd k2(14,1);
    Eigen::ArrayXXd k3(14,1);
    Eigen::ArrayXXd k4(14,1);

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
    Matrice.conservativeResize(Matrice.rows(),Matrice.cols()+2);
    Matrice.col(Matrice.cols()-2)=th.transpose();
    Matrice.col(Matrice.cols()-1)=time.transpose();
    std::ofstream theFile;
    theFile.open(name + "_file.csv" , std::ios::out | std::ios::app);
    theFile << "Position_1,Position_2,Position_3,Velocity_1,Velocity_2,Velocity_3,Mass,APosition_1,APosition_2,APosition_3,AVelocity_1,AVelocity_2,AVelocity_3,AMass,Thrust,Time" << std::endl;
    theFile << Matrice.format(csv)<< std::endl<< std::endl<< std::endl;
    theFile.close();
}

void Indirect_BVP_GA::save(std::string name){
    propagate_save(x0, t0, tGA1, name);
    propagate_save(xGA1_f, tGA1, tGA2, name);
    propagate_save(xGA2_f, tGA2, tf, name);
}
void Indirect_BVP_GA::print(std::string value){
    if (value.compare("Vinf")==0)
    {
        std::cout<<"Vinf Magnitude"<<v_inf_0<<std::endl;
        std::cout<<"Vinf Vector"<<V_inf_0.transpose()/V_inf_0.norm()<<std::endl;
    }
    else if (value.compare("Vot0")==0)
    {
        std::cout<<"Vot0 Vector"<<VOt0.transpose()/VOt0.norm()<<std::endl;
    }
    else if (value.compare("all scalars")==0){
        std::cout<<"m0"<<m0<<std::endl;
        std::cout<<"mLEO"<<mLEO<<std::endl;
        std::cout<<"ve"<<ve<<std::endl;
        std::cout<<"vc"<<vc<<std::endl;
        std::cout<<"Delta V"<<del_v<<std::endl;
    }
}