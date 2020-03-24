//
//  ReadData.cpp
//  
//
//  Created by Clare Cheng on 3/6/20.
//

#include "ReadData.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include "MembLJ.h"
#include <iomanip>


int main(){
    
    int N = 0;
    int NP = 6;
    int Ntot= (N+1)*(N+1);
    Vector in(Ntot+3*NP);
    Vector para(3);
    Vector discPara(2);
    double y;
    
    discPara(0) = N; //Lmax
    discPara(1) = NP; //Number of particles
    

    ifstream myfile;
    myfile.open("Data.txt");//, ios_base::in
    while (!myfile.eof()){
        
        for (int i=0; i<Ntot; i++){
            myfile >> y;
            in(i) = y;
        }
        for (int i=0; i< NP*3; i++){
            myfile >> y;
            in(Ntot + i) = y;
        }
        
        for (int i=0; i<3; i++){
            myfile >> y;
            para(i) = y;
        }
        
         // Construct the model
        MembLJ prob(in, para, discPara);
        prob.f();
        ofstream outFile("Output.txt", ofstream::app);
        outFile.setf(ios_base::scientific);
        //outFile << setprecision(9);
        outFile << para(0) << " \t " << para(1) << " \t "<< para(2) << " \t ";
        outFile << prob._f << endl;
        outFile.close();
        
        
        
    }
    myfile.close();
    
    // Open an output file "output.txt"
    
    // while loop: read each line in Data.txt
    

    

    
    // print to output.txt : parameters and energy (as columns)
    
    // end while loop
    
    
    
        
        
    
    
}
