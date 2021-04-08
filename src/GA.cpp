#include <iostream>
#include <cmath>
#include <queue>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <iomanip>
#include "GA.h"
#include "indirect_method.h"
#include "orbit.h"

Eigen::Vector3i top3(Eigen::VectorXd a){
    int k=3;
    Eigen::Vector3i output;
    std::priority_queue<std::pair<double, int>> q;
    for (int i = 0; i < a.size(); ++i) {
        q.push(std::pair<double, int>(a(i), i));
    }
    for (int i = 0; i < k; ++i) {
        int ki = q.top().second;
        output(i)= ki;
        q.pop();
    }
    return output;
}

//random number generator
double fRand(double fMin, double fMax){
    //srand(time(0));   %was making search space go wrong
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//Crossover Probabity
double beta_dash(double u){
    if (u<=0.5){
        return std::pow(2*u,1.0/3.0);
    }
    else{
        return std::pow(2-2*u, -1.0/3.0);
    }
}

GA::GA(){
    
    Eigen::ArrayXXd FirstGen_param(8,pop_size);
    Eigen::VectorXd FirstGen_fitness(pop_size);
    
    Eigen::VectorXd UpperBound(8);
    Eigen::VectorXd LowerBound(8);

    UpperBound<< t0_up, m0_up, vinf_up, vinf_up, lambdaR0_ub, lambdaR0_ub, lambdaR0_ub, tf_up;
    LowerBound<< t0_lb, m0_lb, vinf_lp, vinf_lp, lambdaR0_lb, lambdaR0_lb, lambdaR0_lb, tf_lb;

    //intialise first generation
    for (int i = 0; i < pop_size; i++){
        FirstGen_param(0,i)=fRand(t0_lb,t0_up);
        FirstGen_param(1,i)=fRand(m0_lb,m0_up);
        FirstGen_param(2,i)=fRand(vinf_lp,vinf_up);
        FirstGen_param(3,i)=fRand(vinf_lp,vinf_up);
        FirstGen_param(4,i)=fRand(lambdaR0_lb,lambdaR0_ub);
        FirstGen_param(5,i)=fRand(lambdaR0_lb,lambdaR0_ub);
        FirstGen_param(6,i)=fRand(lambdaR0_lb,lambdaR0_ub);
        FirstGen_param(7,i)=fRand(tf_lb,tf_up);
        Eigen::VectorXd Pop_param;
        Pop_param=FirstGen_param.col(i);
        IndirectMethod DM(Pop_param);
        FirstGen_fitness(i)=DM.getFitness();
    }
    //GA operators for each generation
    Eigen::ArrayXXd ithGen_param(8,pop_size);
    Eigen::VectorXd ithGen_fitness(pop_size);

    Eigen::ArrayXXd ParentGen_param(8,pop_size);
    Eigen::VectorXd ParentGen_fitness(pop_size);

    Eigen::ArrayXXd i1thGen_param(8,pop_size);
    Eigen::VectorXd i1thGen_fitness(pop_size);

    //array for keeping track of pairings
    std::vector<int> pairings;
    for (int i = 0; i < pop_size; ++i){
        pairings.push_back(i);
    }
    unsigned seed=time(0);

    //set ith Gen as First Gen
    ithGen_param=FirstGen_param;
    ithGen_fitness=FirstGen_fitness;

    for (int i = 0; i < noi; i++){
        //elitisim
        Eigen::Vector3i indices_top3=top3(-ithGen_fitness);
        
        //Tournament Selection
        //1st Tournament
        std::shuffle(pairings.begin(),pairings.end(),std::default_random_engine(seed));
        for (int k = 0; k < pop_size/2; k++){
            if (ithGen_fitness(pairings[k])<ithGen_fitness(pairings[k+pop_size/2])){
                ParentGen_param.col(k)=ithGen_param.col(pairings[k]);
                ParentGen_fitness(k)=ithGen_fitness(pairings[k]);
            }
            else{
                ParentGen_param.col(k)=ithGen_param.col(pairings[k+pop_size/2]);
                ParentGen_fitness(k)=ithGen_fitness(pairings[k+pop_size/2]);
            }    
        }
        //2nd Tournament
        std::shuffle(pairings.begin(),pairings.end(),std::default_random_engine(seed));
        for (int k = 0; k < pop_size/2; k++){
            if (ithGen_fitness(pairings[k])<ithGen_fitness(pairings[k+pop_size/2])){
                ParentGen_param.col(pop_size/2+k)=ithGen_param.col(pairings[k]);
                ParentGen_fitness(pop_size/2+k)=ithGen_fitness(pairings[k]);
            }
            else{
                ParentGen_param.col(pop_size/2+k)=ithGen_param.col(pairings[k+pop_size/2]);
                ParentGen_fitness(pop_size/2+k)=ithGen_fitness(pairings[k+pop_size/2]);
            }    
        }
        std::cout<<i+1<<" out of "<<noi<<" selection done"<<std::endl;

        //Crossover
        std::shuffle(pairings.begin(),pairings.end(),std::default_random_engine(seed));
        for (int k = 0; k < pop_size/2; k++){
            for (int j = 0; j < 8; j++){
                double u=fRand(0,1);
                double beta;
                beta=beta_dash(u);
                double x1, x2;
                x1=ParentGen_param(j,pairings[k]);
                x2=ParentGen_param(j,pairings[k+pop_size/2]);
                i1thGen_param(j,k)=0.5*(x1 + x2 - beta*abs(x2-x1));
                i1thGen_param(j,k+pop_size/2)=0.5*(x1 + x2 + beta*abs(x2-x1));
                // Get the parameters back between the bounds
                if (i1thGen_param(j,k)<LowerBound(j)){
                    i1thGen_param(j,k)=LowerBound(j);
                }
                else if (i1thGen_param(j,k)>UpperBound(j)){
                    i1thGen_param(j,k)=UpperBound(j);
                }
                else if (i1thGen_param(j,k+pop_size/2)<LowerBound(j)){
                    i1thGen_param(j,k+pop_size/2)=LowerBound(j);
                }
                else if (i1thGen_param(j,k+pop_size/2)>UpperBound(j)){
                    i1thGen_param(j,k+pop_size/2)=UpperBound(j);
                }
            }
            Eigen::VectorXd Pop_param1, Pop_param2;
            Pop_param1=i1thGen_param.col(k);
            Pop_param2=i1thGen_param.col(k+pop_size/2);
            IndirectMethod DM1(Pop_param1);
            IndirectMethod DM2(Pop_param2);
            i1thGen_fitness(k)=DM1.getFitness();
            i1thGen_fitness(k+pop_size/2)=DM2.getFitness();
        } 
        std::cout<<i+1<<" out of "<<noi<<" crossover done"<<std::endl;

        //Mutation
        for (int k = 0; k < pop_size; k++){
            double mut_prob=fRand(0,100);
            int check=0;
            if(mut_prob<1){
                i1thGen_param(0,k)=fRand(t0_lb,t0_up);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(1,k)=fRand(m0_lb,m0_up);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(2,k)=fRand(vinf_lp,vinf_up);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(3,k)=fRand(vinf_lp,vinf_up);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(4,k)=fRand(lambdaR0_lb,lambdaR0_ub);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(5,k)=fRand(lambdaR0_lb,lambdaR0_ub);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(6,k)=fRand(lambdaR0_lb,lambdaR0_ub);
                check++;
            }
            mut_prob=fRand(0,100);
            if(mut_prob<1){
                i1thGen_param(7,k)=fRand(tf_lb,tf_up);
                check++;
            }
            if (check!=0){
                Eigen::VectorXd Pop_param3;
                Pop_param3=i1thGen_param.col(k);
                IndirectMethod DM3(Pop_param3);
                i1thGen_fitness(k)=DM3.getFitness();
            }
        }
        std::cout<<i+1<<" out of "<<noi<<" mutation done"<<std::endl;

        //elitism part 2
        //Finding worst 3 fitness indices for i+1 th Gen
        Eigen::Vector3i indices_bot3=top3(i1thGen_fitness);
        for (int k = 0; k < 3; k++){
            i1thGen_fitness(indices_bot3(k))=ithGen_fitness(indices_top3(k));
            i1thGen_param.col(indices_bot3(k))=ithGen_param.col(indices_top3(k));
        }
        ithGen_param=i1thGen_param;
        ithGen_fitness=i1thGen_fitness;
        std::cout<<i+1<<" out of iterations "<<noi<<" done"<<std::endl;
    }
    BestGene=Eigen::VectorXd::Zero(8);
    Eigen::VectorXd::Index min_index;
    BestFitness=ithGen_fitness.minCoeff(&min_index);
    BestGene=ithGen_param.col(min_index);
}

double GA::get_BestFitness(){
    return BestFitness;
}
Eigen::VectorXd GA::get_BestGene(){
    return BestGene;
}