#include <iostream>
#include <cmath>
#include <queue>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <iomanip>
#include "Indirect_BVP_GA.h"
#include "Genetic_GA.h"
#include "orbit.h"

Genetic_GA::Genetic_GA(std::string seq){
    
    Eigen::ArrayXXd FirstGen_param(24,pop_size);
    Eigen::VectorXd FirstGen_fitness(pop_size);
    
    Eigen::VectorXd UpperBound(24);
    Eigen::VectorXd LowerBound(24);
    LowerBound(0)=t0_lb;
    LowerBound(1)=m0_lb;
    LowerBound(2)=vinf_dep_lp;
    LowerBound(3)=vinf_dep_lp;
    LowerBound(4)=lambdaR0_lb;
    LowerBound(5)=lambdaR0_lb;
    LowerBound(6)=lambdaR0_lb;
    // First Gravity Assist Parameters
    LowerBound(7)=tGA1_lb;
    LowerBound(8)=vinf_GA_lp;
    LowerBound(9)=vinf_GA_lp;
    LowerBound(10)=lambdaR0_lb;
    LowerBound(11)=lambdaR0_lb;
    LowerBound(12)=lambdaR0_lb;
    // Second Gravity Assist Parameters (with minimum height constraint)
    LowerBound(13)=tGA2_lb;
    LowerBound(14)=vinf_GA_lp;
    LowerBound(15)=vinf_GA_lp;
    LowerBound(16)=lambdaR0_lb;
    LowerBound(17)=lambdaR0_lb;
    LowerBound(18)=lambdaR0_lb;
    LowerBound(19)=lambdaR0_lb;
    LowerBound(20)=lambdaR0_lb;
    LowerBound(21)=lambdaR0_lb;
    // Destination/ Final Parameters
    LowerBound(22)=tf_lb;
    LowerBound(23)=tstar_lb;

    UpperBound(0)=t0_up;
    UpperBound(1)=m0_up;
    UpperBound(2)=vinf_dep_up;
    UpperBound(3)=vinf_dep_up;
    UpperBound(4)=lambdaR0_ub;
    UpperBound(5)=lambdaR0_ub;
    UpperBound(6)=lambdaR0_ub;
    // First Gravity Assist Parameters
    UpperBound(7)=tGA1_up;
    UpperBound(8)=vinf_GA_up;
    UpperBound(9)=vinf_GA_up;
    UpperBound(10)=lambdaR0_ub;
    UpperBound(11)=lambdaR0_ub;
    UpperBound(12)=lambdaR0_ub;
    // Second Gravity Assist Parameters (with minimum height constraint)
    UpperBound(13)=tGA2_up;
    UpperBound(14)=vinf_GA_up;
    UpperBound(15)=vinf_GA_up;
    UpperBound(16)=lambdaR0_ub;
    UpperBound(17)=lambdaR0_ub;
    UpperBound(18)=lambdaR0_ub;
    UpperBound(19)=lambdaR0_ub;
    UpperBound(20)=lambdaR0_ub;
    UpperBound(21)=lambdaR0_ub;
    // Destination/ Final Parameters
    UpperBound(22)=tf_up;
    UpperBound(23)=tstar_up;

    //intialise first generation
    for (int i = 0; i < pop_size; i++){
        for (int j = 0; j < 24; j++){
            FirstGen_param(j,i)=fRand(LowerBound(j),UpperBound(j));
        }
        Eigen::VectorXd Pop_param;
        Pop_param=FirstGen_param.col(i);
        Indirect_BVP_GA GA(Pop_param, seq);
        FirstGen_fitness(i)=GA.getFitness();
    }

    //GA operators for each generation
    Eigen::ArrayXXd ithGen_param(24,pop_size);
    Eigen::VectorXd ithGen_fitness(pop_size);

    Eigen::ArrayXXd ParentGen_param(24,pop_size);
    Eigen::VectorXd ParentGen_fitness(pop_size);

    Eigen::ArrayXXd i1thGen_param(24,pop_size);
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
            for (int j = 0; j < 24; j++){
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
            Indirect_BVP_GA GA1(Pop_param1, seq);
            Indirect_BVP_GA GA2(Pop_param2, seq);
            i1thGen_fitness(k)=GA1.getFitness();
            i1thGen_fitness(k+pop_size/2)=GA2.getFitness();
        } 
        std::cout<<i+1<<" out of "<<noi<<" crossover done"<<std::endl;

        //Mutation
        for (int k = 0; k < pop_size; k++){
            double mut_prob;
            int check=0;
            for (int j = 0; j < 24; j++){
                mut_prob=fRand(0,100);
                if(mut_prob<1){
                    i1thGen_param(j,k)=fRand(LowerBound(j),UpperBound(j));
                    check++;
                }
            }
            if (check!=0){
                Eigen::VectorXd Pop_param3;
                Pop_param3=i1thGen_param.col(k);
                Indirect_BVP_GA GA3(Pop_param3,seq);
                i1thGen_fitness(k)=GA3.getFitness();
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
    BestGene=Eigen::VectorXd::Zero(24);
    Eigen::VectorXd::Index min_index;
    BestFitness=ithGen_fitness.minCoeff(&min_index);
    BestGene=ithGen_param.col(min_index);
}

double Genetic_GA::get_BestFitness(){
    return BestFitness;
}
Eigen::VectorXd Genetic_GA::get_BestGene(){
    return BestGene;
}

Eigen::Vector3i Genetic_GA::top3(Eigen::VectorXd a){
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
double Genetic_GA::fRand(double fMin, double fMax){
    //srand(time(0));   %was making search space go wrong
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//Crossover Probabity
double Genetic_GA::beta_dash(double u){
    if (u<=0.5){
        return std::pow(2*u,1.0/3.0);
    }
    else{
        return std::pow(2-2*u, -1.0/3.0);
    }
}