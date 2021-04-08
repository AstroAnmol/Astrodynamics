#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <eigen-3.3.7/Eigen/Dense>
#include <cmath>
#include <queue>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include "orbit.h"
#include "cr3bp.h"
#include "indirect_method.h"
#include "GA.h"


int main(){
        srand(time(0));
        /*
        Eigen::VectorXd Gene(8);
        Gene<<50.6555, 0.972694, -0.000650264, -0.000551347, -0.590737, -1.39754, -0.552314, 63.0691;
        IndirectMethod Trial(Gene);
        Trial.print("VEt0");
        Trial.print("Vinf");
        Trial.print("all scalars");
        */
        
        for (int i = 0; i < 10; i++)
        {
            GA trial1;
            Eigen::VectorXd BestGeneT1;
            BestGeneT1=trial1.get_BestGene();;
            std::cout<<BestGeneT1.transpose()<<std::endl;
            std::cout<<trial1.get_BestFitness()<<std::endl;
            IndirectMethod Best(BestGeneT1);
            std::string name;
            std::ostringstream oss;
            oss << "Trial_" << i;
            name=oss.str();
            Best.save(name);
        }
        
}