import numpy as np

def differential(x): #x is the following vector: R, V, m, Lambda_R, Lambda_V, lambda_m
    R=x[:3]
    V=x[3:6]
    m=x[6]
    Lambda_R=x[7:10]
    Lambda_V=x[10:13]
    lambda_m=x[13]

    r=np.linalg.norm(R)
    lambda_V=np.linalg.norm(Lambda_V)

    T=Thrust(x)
    t=np.linalg.norm(T)

    output=np.zeros(14)
    output[0]=V[0]
    output[1]=V[1]
    output[2]=V[2]
    output[3]=-R[0]/r**3 + T[0]/m
    output[4]=-R[1]/r**3 + T[1]/m
    output[5]=-R[2]/r**3 + T[2]/m
    output[6]=-t
    output[7]=-Lambda_V[0]/r**3 + 3*R[0]*(np.dot(R,Lambda_V))/r**5
    output[8]=-Lambda_V[1]/r**3 + 3*R[1]*(np.dot(R,Lambda_V))/r**5
    output[9]=-Lambda_V[2]/r**3 + 3*R[2]*(np.dot(R,Lambda_V))/r**5
    output[10]=-Lambda_R[0]
    output[11]=-Lambda_R[1]
    output[12]=-Lambda_R[2]
    output[13]=t*lambda_V/m**2

    return output

def Thrust(x):
    R=x[:3]
    V=x[3:6]
    m=x[6]
    Lambda_R=x[7:10]
    Lambda_V=x[10:13]
    lambda_m=x[13]

    r=np.linalg.norm(R)
    lambda_V=np.linalg.norm(Lambda_V)

    T0=0.02
    Tmax=T0/r**2
    
    SF=lambda_V/m - lambda_m

    if SF>0:
        return Tmax*Lambda_V/lambda_V
    elif SF<0:
        return np.zeros(3)

def Kepler_prob(R0, V0,del_t):
    alpha=1/a
    chi_0
    chi=chi_0
    tol=1e-8
    ratio=1
    r0=np.linalg.norm(R0)
    v0=np.linalg.norm(V0)
    if e>1:
        x=-2*mu*del_t/(a*(R0.dot(V0) + del_t * np.sqrt(-a*mu) * (1-r0/a)))
        chi_0=del_t*std::sqrt(-a)*std::log(x)
    }
    else {
        chi_0= std::sqrt(mu) * std::abs(alpha) *del_t;
    }
    int i=0;
    int maxiter=10000;
    while(std::abs(ratio)>tol){
        if(i>maxiter){
            break;}
        z=chi*chi*alpha;
        double f_chi, f_dash_chi;
        f_chi=R0.dot(V0)*chi*chi*stumpff_C(z)/(std::sqrt(mu));
        f_chi=f_chi+r0*chi*(1-z*stumpff_S(z));
        f_chi=f_chi+ std::pow(chi,3)*stumpff_S(z)- std::sqrt(mu)*del_t;

        f_dash_chi=R0.dot(V0)*chi*(1-z*stumpff_S(z))/(std::sqrt(mu));
        f_dash_chi=f_dash_chi + chi*chi*stumpff_C(z) + r0 -z*r0*stumpff_C(z);

        ratio=f_chi/f_dash_chi;
        chi=chi-ratio;
        i++;
    }
    
    z=chi*chi*alpha;
    double f, g, f_dot, g_dot;
    Eigen::Vector3d R_t, V_t;
    f= 1- chi*chi*stumpff_C(z)/r0;
    g=del_t - (chi*chi*chi*stumpff_S(z))/std::sqrt(mu);

    R_t=f*R0 +g*V0;
    double r_t=R_t.norm();
    f_dot=(std::sqrt(mu)/(r0*r_t)) * chi*(z*stumpff_S(z)-1);
    g_dot=1- chi*chi*stumpff_C(z)/r_t;

    V_t=f_dot*R0 + g_dot*V0;
    Eigen::VectorXd res(6);
    res.block(0,0,3,1)=R_t;
    res.block(3,0,3,1)=V_t;
    return res;
 }