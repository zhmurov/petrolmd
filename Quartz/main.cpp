#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "xyzio.h"
#include "psfio.h"
#include "pdbio.h"

#define NX 2
#define NY 2

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
    while (std::getline(rtpInStream, line) && line.compare("[ SiO ]") != 0)
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

    PDB pdbOut;
    pdbOut.atomCount = atomsOut.size();
    pdbOut.atoms = (PDBAtom*)calloc(pdbOut.atomCount, sizeof(PDBAtom));
    pdbOut.ssCount = 0;
    pdbOut.symmetryCount = 0;
    pdbOut.matrixCount = 0;
    idx = 0;
    for (auto atomOut : atomsOut)
    {
        pdbOut.atoms[idx].id = idx + 1;
        sprintf(pdbOut.atoms[idx].name, "%s", atomOut.name.c_str());
        pdbOut.atoms[idx].chain = 'A';
        sprintf(pdbOut.atoms[idx].resName, "SiO");
        pdbOut.atoms[idx].altLoc = ' ';
        pdbOut.atoms[idx].resid = atomOut.ix + NX*atomOut.iy;
        pdbOut.atoms[idx].x = atomOut.x;
        pdbOut.atoms[idx].y = atomOut.y;
        pdbOut.atoms[idx].z = atomOut.z;
        sprintf(pdbOut.atoms[idx].segment, "A");
        pdbOut.atoms[idx].occupancy = 0.0;
        pdbOut.atoms[idx].beta = 0.0;
        idx ++;
    }
    writePDB("slab.pdb", &pdbOut);
}