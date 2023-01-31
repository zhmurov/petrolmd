#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "xyzio.h"

#define NX 10
#define NY 10
#define LZ 200.0

#define DOPBC true

struct Atom
{
    std::string name;
    std::string type;
    float charge;
    int cgr;
    float x, y, z;
    int ix;
    int iy;
};

int main(int argc, char *argv[])
{
    XYZ xyzIn;
    readXYZ("/home/zhmurov/git/artemzhmurov/petrolmd/Quartz/files/input.xyz", &xyzIn);
    
    std::ifstream ucFile;
    ucFile.open("/home/zhmurov/git/artemzhmurov/petrolmd/Quartz/files/crystal.dat");

    float a, b, c, alpha, beta, gamma;

    ucFile >> a;
    ucFile >> b;
    ucFile >> c;
    ucFile >> alpha;
    ucFile >> beta;
    ucFile >> gamma;

    std::cout << a << "  " << b << "  " << c << "  " << alpha << "  " << beta << "  " << gamma << std::endl;

    std::cout << "Recommended x - y box size is " << a*NX << " - " << b*NY << std::endl;

    ucFile.close();

    XYZ xyzOut;
    xyzOut.atomCount = xyzIn.atomCount*NX*NY;
    xyzOut.atoms = (XYZAtom*)calloc(xyzOut.atomCount, sizeof(XYZAtom));
    for (int iy = 0; iy < NY; iy ++)
    {
        for (int ix = 0; ix < NX; ix ++)
        {
            for (int i = 0; i < xyzIn.atomCount; i++)
            {
                int iOut = i + ix*xyzIn.atomCount + iy*NX*xyzIn.atomCount;
                xyzOut.atoms[iOut].name = xyzIn.atoms[i].name;
                float x = xyzIn.atoms[i].x;
                float y = xyzIn.atoms[i].y;
                float z = xyzIn.atoms[i].z;

                xyzOut.atoms[iOut].x = x + a*ix;
                xyzOut.atoms[iOut].y = y + b*iy;
                xyzOut.atoms[iOut].z = z;

            }
        }
    }
    writeXYZ("slab.xyz", &xyzOut);

    std::vector<Atom> atomsIn;

    std::fstream rtpInStream;
    rtpInStream.open("/home/zhmurov/git/artemzhmurov/charmm36/charmm36.ff/silicates.rtp");

    std::string line;
    while (std::getline(rtpInStream, line) && line.compare("[ Q011 ]") != 0)
    {
        
    }
    while (std::getline(rtpInStream, line) && line.compare("  [ atoms ]") != 0)
    {
    }
    while (std::getline(rtpInStream, line) && line.compare("  [ bonds ]") != 0)
    {
        //std::cout << line << std::endl;
        std::stringstream lineStream(line); 
        Atom atom;
        lineStream >> atom.name >> atom.type >> atom.charge >> atom.cgr;
        atomsIn.push_back(atom);

    }
    rtpInStream.close();

    std::cout << "Atoms:" << std::endl;
    for (auto atom : atomsIn)
    {
        std::cout << atom.name << " " << atom.type << " " << atom.charge << " " << atom.cgr << std::endl;
    }
    
    std::vector<Atom> atomsOut;

    int idx = 0;

    for (int iy = 0; iy < NY; iy ++)
    {
        for (int ix = 0; ix < NX; ix ++)
        {
            for (auto atomIn : atomsIn)
            {
                Atom atomOut;
                atomOut.name = atomIn.name;
                atomOut.type = atomIn.type;
                atomOut.charge = atomIn.charge;
                atomOut.cgr = atomIn.cgr;
                atomOut.ix = ix;
                atomOut.iy = iy;
                atomsOut.push_back(atomOut);
                idx++;
            }
        }
    }

    for (int i = 0; i < xyzOut.atomCount; i++)
    {
        atomsOut[i].x = xyzOut.atoms[i].x;
        atomsOut[i].y = xyzOut.atoms[i].y;
        atomsOut[i].z = xyzOut.atoms[i].z;
    }

    FILE* groOut = fopen("slab.gro", "w");
    fprintf(groOut, "Quartz slab (%dx%d)\n", NX, NY);
    fprintf(groOut, "%d\n", static_cast<int>(atomsOut.size()));
    idx = 0;
    for (auto atomOut : atomsOut)
    {
        if (DOPBC)
        {
            atomOut.x -= 0.5*a;
            atomOut.y -= 0.5*b;
            if (atomOut.x < 0)
            {
                atomOut.x += a*NX;
            }
            if (atomOut.y < 0)
            {
                atomOut.y += b*NY;
            }
            atomOut.z += 0.5*LZ;
        }
        fprintf(groOut, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
            atomOut.ix + NX*atomOut.iy + 1,
            "Q011",
            atomOut.name.c_str(),
            idx + 1,
            atomOut.x*0.1, // A to nm
            atomOut.y*0.1,
            atomOut.z*0.1);
        idx++;
    }
    fprintf(groOut, "%f %f %f\n", a*NX*0.1, b*NY*0.1, LZ*0.1);
    fclose(groOut);
}