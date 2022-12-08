#include <fstream>
#include <iostream>
#include <string>
#include <iostream>

#include "xyzio.h"

#define N_X 2
#define N_Y 2 

int main(int argc, char *argv[])
{
    XYZ xyzIn;
    readXYZ("Quartz/files/input.xyz", &xyzIn);
    
    std::ifstream ucFile;
    ucFile.open("Quartz/files/crystal.dat");

    float a, b, c, alpha, beta, gamma;

    ucFile >> a;
    ucFile >> b;
    ucFile >> c;
    ucFile >> alpha;
    ucFile >> beta;
    ucFile >> gamma;

    std::cout << a << "  " << b << "  " << c << "  " << alpha << "  " << beta << "  " << gamma << std::endl;

    ucFile.close();


    XYZ xyzOut;
    xyzOut.atomCount = xyzIn.atomCount*N_X*N_Y;
}