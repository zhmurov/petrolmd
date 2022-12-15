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

struct Bond
{
    std::string name1;
    std::string name2;
};

std::map<std::string, int> atomNameToIdx;

std::string nameToShortName(std::string name)
{
    std::stringstream stream;
    stream << std::hex << std::setfill('0') << std::setw(3) << atomNameToIdx[name];
    return stream.str();
}

std::string convertName(std::string name, int ix, int iy)
{
    return name + "_" + std::to_string(ix) + "_" + std::to_string(iy);
}

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
    std::vector<Bond> bondsIn;

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
    while (std::getline(rtpInStream, line) && !line.empty())
    {
        //std::cout << line << std::endl;
        std::stringstream lineStream(line); 
        Bond bond;
        lineStream >> bond.name1 >> bond.name2;
        bondsIn.push_back(bond);
    }
    rtpInStream.close();

    std::cout << "Atoms:" << std::endl;
    for (auto atom : atomsIn)
    {
        std::cout << atom.name << " " << atom.type << " " << atom.charge << " " << atom.cgr << std::endl;
    }
    
    std::cout << "Bonds:" << std::endl;
    for (auto bond : bondsIn)
    {
        std::cout << bond.name1 << " " << bond.name2 << std::endl;
    }

    std::vector<Atom> atomsOut;
    std::vector<Bond> bondsOut;

    int idx = 0;

    for (int iy = 0; iy < NY; iy ++)
    {
        for (int ix = 0; ix < NX; ix ++)
        {
            for (auto atomIn : atomsIn)
            {
                Atom atomOut;
                atomOut.name = convertName(atomIn.name, ix, iy);
                atomOut.type = atomIn.type;
                atomOut.charge = atomIn.charge;
                atomOut.cgr = atomIn.cgr;
                atomOut.ix = ix;
                atomOut.iy = iy;
                atomsOut.push_back(atomOut);
                atomNameToIdx.insert({atomOut.name, idx});
                idx++;
            }
        }
    }
    for (int iy = 0; iy < NY; iy ++)
    {
        for (int ix = 0; ix < NX; ix ++)
        {
            for (auto bondIn : bondsIn)
            {
                Bond bondOut;
                bondOut.name1 = convertName(bondIn.name1, ix, iy);
                bondOut.name2 = convertName(bondIn.name2, ix, iy);
                bondsOut.push_back(bondOut);
            }
            //if (ix > 0)
            {
                int ixm1 = (ix > 0) ? ix-1 : NX-1;
                Bond bond1;
                bond1.name1 = convertName("SI4", ixm1, iy);
                bond1.name2 = convertName("O12", ix, iy);
                bondsOut.push_back(bond1);

                Bond bond2;
                bond2.name1 = convertName("SI27", ixm1, iy);
                bond2.name2 = convertName("O35", ix, iy);
                bondsOut.push_back(bond2);

                Bond bond3;
                bond3.name1 = convertName("SI22", ixm1, iy);
                bond3.name2 = convertName("O18", ix, iy);
                bondsOut.push_back(bond3);

                Bond bond4;
                bond4.name1 = convertName("O9", ixm1, iy);
                bond4.name2 = convertName("SI25", ix, iy);
                bondsOut.push_back(bond4);
            }
            //if (iy > 0)
            {
                int iym1 = (iy > 0) ? iy-1 : NY-1 ;
                Bond bond1;
                bond1.name1 = convertName("SI1", ix, iym1);
                bond1.name2 = convertName("O17", ix, iy);
                bondsOut.push_back(bond1);

                Bond bond2;
                bond2.name1 = convertName("O30", ix, iym1);
                bond2.name2 = convertName("SI27", ix, iy);
                bondsOut.push_back(bond2);

                Bond bond3;
                bond3.name1 = convertName("O31", ix, iym1);
                bond3.name2 = convertName("SI28", ix, iy);
                bondsOut.push_back(bond3);

                Bond bond4;
                bond4.name1 = convertName("SI23", ix, iym1);
                bond4.name2 = convertName("O7", ix, iy);
                bondsOut.push_back(bond4);
            }
        }
    }

    for (int i = 0; i < xyzOut.atomCount; i++)
    {
        atomsOut[i].x = xyzOut.atoms[i].x;
        atomsOut[i].y = xyzOut.atoms[i].y;
        atomsOut[i].z = xyzOut.atoms[i].z;
    }

    std::ofstream rtpOutFile;
    rtpOutFile.open("slab.rtp");

    rtpOutFile << "[ bondedtypes ]" << std::endl;
    rtpOutFile << "; bonds  angles  dihedrals  impropers  all_dihedrals  nrexcl  HH14  RemoveDih " << std::endl;
    rtpOutFile << "1       5        9          2            1           3      1       0" << std::endl;

    //rtpOutFile << "[ SiO_" << NX << "_" << NY << " ]" << std::endl;
    rtpOutFile << "[ SiO ]" << std::endl;
    rtpOutFile << "; " << std::endl;
    rtpOutFile << "[ atoms ]" << std::endl;
    for (auto atomOut : atomsOut)
    {
        //rtpOutFile << "  " << atomOut.name << "     " << atomOut.type << " " << atomOut.charge << "   " << atomOut.cgr << std::endl;
        rtpOutFile << "  " << nameToShortName(atomOut.name) << "     " << atomOut.type << " " << atomOut.charge << "   " << atomOut.cgr << std::endl;
    }
    rtpOutFile << "[ bonds ]" << std::endl;
    for (auto bondOut : bondsOut)
    {
        //rtpOutFile << "   " << bondOut.name1 << "   " << bondOut.name2 << std::endl;
        rtpOutFile << "   " << nameToShortName(bondOut.name1) << "   " << nameToShortName(bondOut.name2) << std::endl;
    }
    rtpOutFile.close();

    PSF psfOut;
    psfOut.natom = atomsOut.size();
    psfOut.nbond = bondsOut.size();
    psfOut.ntheta = 0;
    psfOut.nphi = 0;
    psfOut.nimphi = 0;
    psfOut.ncmap = 0;
    psfOut.nnb = 0;

    psfOut.atoms = (PSFAtom*)calloc(psfOut.natom, sizeof(PSFAtom));
    psfOut.bonds = (PSFBond*)calloc(psfOut.nbond, sizeof(PSFBond));

    idx = 0;
    for (auto atomOut : atomsOut)
    {
        PSFAtom atom;
        atom.id = atomNameToIdx[atomOut.name] + 1;
        sprintf(atom.name, "%3s ", atomOut.name.substr(0, 3).c_str());
        sprintf(atom.type, "%s", atomOut.type.c_str());
        atom.q = atomOut.charge;
        
        atom.resid = 1 + atomOut.ix + NX*atomOut.iy;
        sprintf(atom.resName, "SiO2");
        sprintf(atom.segment, "SiO2");
        atom.m = 1.0;
        psfOut.atoms[idx] = atom;
        idx ++;
    }
    idx = 0;
    for (auto bondOut : bondsOut)
    {
        PSFBond bond;
        bond.i = atomNameToIdx[bondOut.name1] + 1;
        bond.j = atomNameToIdx[bondOut.name2] + 1;
        psfOut.bonds[idx] = bond;
        idx ++;
    }
    writePSF("slab.psf", &psfOut);

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
        sprintf(pdbOut.atoms[idx].name, "%s", nameToShortName(atomOut.name).c_str());
        pdbOut.atoms[idx].chain = 'A';
        sprintf(pdbOut.atoms[idx].resName, "SiO");
        pdbOut.atoms[idx].altLoc = ' ';
        pdbOut.atoms[idx].resid = 1;
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