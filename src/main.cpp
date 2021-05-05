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
#include "Indirect_BVP_DM.h"
#include "Genetic_DM.h"


int main(){
        srand(time(0));
        /*
        Eigen::VectorXd Gene(8);
        Gene<<50.6555, 0.972694, -0.000650264, -0.000551347, -0.590737, -1.39754, -0.552314, 63.0691;
        IndirectMethod Trial(Gene);
        Trial.print("VEt0");
        Trial.print("Vinf");
        Trial.print("all scalars");
        Trial.save("Try");
        */
        
        for (int i = 0; i < 1; i++)
        {
            Genetic_DM trial1;
            Eigen::VectorXd BestGeneT1;
            BestGeneT1=trial1.get_BestGene();;
            std::cout<<BestGeneT1.transpose()<<std::endl;
            std::cout<<trial1.get_BestFitness()<<std::endl;
            Indirect_BVP_DM Best(BestGeneT1);
            std::string name;
            std::ostringstream oss;
            oss << "DirectMars_" << i;
            name=oss.str();
            Best.save(name);
        }      
}