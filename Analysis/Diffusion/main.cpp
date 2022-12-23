#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <ranges>
#include <chemfiles.hpp>

#include "pdbio.h"

#define NBINX 40
#define NBINY 40
#define NBINZ 40

inline int getIdx(int ix, int iy, int iz)
{
    assert(ix < NBINX);
    assert(iy < NBINY);
    assert(iz < NBINZ);
    assert(ix >= 0);
    assert(iy >= 0);
    assert(iz >= 0);
    return ix + iy*NBINX + iz*NBINX*NBINY;
}

inline int getIdx(int* ixyz)
{
    return getIdx(ixyz[0], ixyz[1], ixyz[2]);
}

std::map<std::string, int> atomCounts = {
    {"C1H4", 5},
    {"C2H6", 8},
    {"C3H8", 11},
    {"C4H10", 14},
    {"C4H10_ISO", 14},
    {"C5H12", 17},
    {"C5H12_ISO", 17},
    {"C6H14", 20},
    {"C7H16", 23},
    {"C8H18", 26},
    {"C9H20", 29},
    {"C10H22", 32},
    {"C11H24", 35},
    {"C12H26", 38},
    {"C13H28", 41},
    {"C14H30", 44},
    {"C15H32", 47},
    {"C16H34", 50},
    {"C17H36", 53},
    {"C18H38", 56},
    {"C19H40", 59},
    {"C20H42", 62},
    {"C21H44", 65},
    {"C22H46", 68},
    {"C23H48", 71},
    {"C24H50", 74},
    {"C25H52", 77},
    {"C26H54", 80},
    {"C27H56", 83},
    {"C28H58", 86},
    {"C29H60", 89},
    {"C30H62", 92},
    {"C31H64", 95},
    {"C32H66", 98},
    {"C33H68", 101},
    {"C34H70", 104},
    {"C35H72", 107},
    {"C36H74", 110},
    {"SOL", 3},
    {"H2O", 3},
    {"NA", 1},
    {"CL", 1}
};

struct MoleculeType
{
    std::string name;
    int atomCount = 0;
    int count;
};

std::vector<MoleculeType> moleculeTypes;
int totMoleculesCount;

void readMoleculeTypes(std::string filename)
{
    std::fstream stream;
    stream.open(filename);

    std::string line;
    while (std::getline(stream, line) && line.compare("[ molecules ]") != 0)
    {
        
    }
    std::getline(stream, line);

    int totAtomCount = 0;
    totMoleculesCount = 0;
    while (std::getline(stream, line) && !line.empty())
    {
        std::cout << line << std::endl;
        std::stringstream lineStream(line); 
        MoleculeType moleculeType;
        lineStream >> moleculeType.name >> moleculeType.count;
        moleculeType.atomCount = atomCounts[moleculeType.name];
        moleculeTypes.push_back(moleculeType);

        totAtomCount += moleculeType.count*moleculeType.atomCount;
        totMoleculesCount += moleculeType.count;
    }
    stream.close();

    std::cout << "Molecules read:" << std::endl;
    for (auto moleculeType : moleculeTypes)
    {
        std::cout << moleculeType.name << " " << moleculeType.count << "  " << moleculeType.atomCount << " atoms each" << std::endl;
    }
    std::cout << "Total number of atoms: " << totAtomCount << std::endl;
}

std::vector<chemfiles::Vector3D> computeCOMs(chemfiles::span<chemfiles::Vector3D> positions)
{
    std::vector<chemfiles::Vector3D> coms = std::vector<chemfiles::Vector3D>(totMoleculesCount, {0.0, 0.0, 0.0});
    int molid = 0;
    int firstAtom = 0;
    for (auto moleculeType : moleculeTypes)
    {
        for (int molnum = 0; molnum < moleculeType.count; molnum++)
        {
            chemfiles::Vector3D com = {0.0, 0.0, 0.0};
            for (int atomid = 0; atomid < moleculeType.atomCount; atomid++)
            {
                com = com + positions[firstAtom + atomid];
            }
            com[0] /= moleculeType.atomCount;
            com[1] /= moleculeType.atomCount;
            com[2] /= moleculeType.atomCount;
            assert(molid < totMoleculesCount);
            coms[molid] = com;
            molid ++;
            firstAtom += moleculeType.atomCount;
        }
    }
    return coms;
}


int main() {
    int nbin[3];
    nbin[0] = NBINX;
    nbin[1] = NBINY;
    nbin[2] = NBINZ;

    readMoleculeTypes("/home/zhmurov/Data/Sirius/PetrolMD/25-charmm-yamburg/yamburg_recomb.top");

    chemfiles::Trajectory mds("/home/zhmurov/Data/Sirius/PetrolMD/25-charmm-yamburg/md_iso.xtc");
    mds.set_topology("/home/zhmurov/Data/Sirius/PetrolMD/25-charmm-yamburg/md_iso.tpr");

    //chemfiles::Trajectory mds("/home/zhmurov/Sirius/5-charmm/md_iso.xtc");
    //mds.set_topology("/home/zhmurov/Sirius/5-charmm/md_iso.tpr");

    std::vector<double> hist = std::vector<double>(NBINX*NBINY*NBINZ, 0.0);

    /*PDB pdbOut;
    pdbOut.atomCount = NBINX*NBINY*NBINZ;
    pdbOut.atoms = (PDBAtom*)calloc(pdbOut.atomCount, sizeof(PDBAtom));
    pdbOut.ssCount = 0;
    pdbOut.symmetryCount = 0;
    pdbOut.matrixCount = 0;
    for (int ix = 0; ix < NBINX; ix++)
    {
        for (int iy = 0; iy < NBINY; iy++)
        {
            for (int iz = 0; iz < NBINZ; iz++)
            {
                int idx = getIdx(ix, iy, iz);
                pdbOut.atoms[idx].id = (idx + 1) % 100000;
                sprintf(pdbOut.atoms[idx].name, "C");
                pdbOut.atoms[idx].chain = 'A';
                sprintf(pdbOut.atoms[idx].resName, "His");
                pdbOut.atoms[idx].altLoc = ' ';
                pdbOut.atoms[idx].resid = 1;
                sprintf(pdbOut.atoms[idx].segment, "A");
                pdbOut.atoms[idx].beta = 0.0;
            }
        }
    }*/
    auto frame = mds.read();
    
    std::vector<chemfiles::Vector3D> oldcoms = std::vector<chemfiles::Vector3D>(totMoleculesCount, {0.0, 0.0, 0.0});
    std::vector<double> dr2s = std::vector<double>(moleculeTypes.size(), 0.0);
    /*std::vector<double> dr2_liquid = std::vector<double>(frame.size(), 0.0);
    std::vector<double> dr2_gas = std::vector<double>(frame.size(), 0.0);
    std::vector<int> ndr2_liquid = std::vector<int>(frame.size(), 0);
    std::vector<int> ndr2_gas = std::vector<int>(frame.size(), 0);*/

    std::cout << "Atom count: " << frame.size() << std::endl;

    for (auto atom : frame.topology())
    {
        
        /*std::cout << atom.name() << ": ";
        for (auto prop : atom.properties().value())
        {
            std::cout << prop.first << "  " << prop.second.as_string();
        }
        std::cout << std::endl;*/

    }

    for (size_t step = 1; step < mds.nsteps(); step++)
    {
        auto frame = mds.read();
        auto uc = frame.cell();
        auto ucl = uc.lengths();
        std::cout << ucl[0] << "  " << ucl[1] << "  " << ucl[2] << std::endl;
        
        int idx[3];
        
        for (auto pos : frame.positions())
        {   
            for (int d = 0; d < 3; d++)
            {
                // TODO: Deal with PBC more elegantly
                if (pos[d] < 0)
                {
                    pos[d] += ucl[d];
                } else
                if (pos[d] >= ucl[d])
                {
                    pos[d] -= ucl[d];
                }
                idx[d] = floor((pos[d]*nbin[d])/ucl[d]);
                //std::cout << idx[d] << "   ";
                assert(idx[d] >= 0);
                assert(idx[d] < nbin[d]);
            }
            //std::cout << std::endl;
            hist[getIdx(idx)] ++;
        }

        /*for (int ix = 0; ix < NBINX; ix++)
        {
            for (int iy = 0; iy < NBINY; iy++)
            {
                for (int iz = 0; iz < NBINZ; iz++)
                {
                    int idx = getIdx(ix, iy, iz);
                    pdbOut.atoms[idx].x = static_cast<float>(ix)*ucl[0]/static_cast<float>(NBINX);
                    pdbOut.atoms[idx].y = static_cast<float>(iy)*ucl[1]/static_cast<float>(NBINY);
                    pdbOut.atoms[idx].z = static_cast<float>(iz)*ucl[2]/static_cast<float>(NBINZ);
                    pdbOut.atoms[idx].occupancy = hist[idx];
                }
            }
        }*/

        std::vector<chemfiles::Vector3D> coms = computeCOMs(frame.positions());

        int molid = 0;
        int moleculeTypeId = 0;
        for (auto moleculeType : moleculeTypes)
        {
            for (int molnum = 0; molnum < moleculeType.count; molnum++)
            {
                chemfiles::Vector3D dr = coms[molid] - oldcoms[molid];
                double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                dr2s[moleculeTypeId] += dr2;
                molid ++;
            }
            moleculeTypeId++;
        }

        std::copy(coms.begin(), coms.end(), std::back_inserter(oldcoms));

        /*int i = 0;
        for (auto pos : frame.positions())
        {   
            for (int d = 0; d < 3; d++)
            {
                // TODO: Deal with PBC more elegantly
                if (pos[d] < 0)
                {
                    pos[d] += ucl[d];
                } else
                if (pos[d] >= ucl[d])
                {
                    pos[d] -= ucl[d];
                }
                idx[d] = floor((pos[d]*nbin[d])/ucl[d]);
                assert(idx[d] >= 0);
                assert(idx[d] < nbin[d]);
            }
            chemfiles::Vector3D dr = pos - oldpos[i];
            double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            if (hist[getIdx(idx)] > 5)
            {
                dr2_liquid[i] += dr2;// in liquid
                ndr2_liquid[i] ++;
            }
            else
            {
                dr2_gas[i] += dr2;// in gas
                ndr2_gas[i]++;
            }
            i++;
        }
        std::copy(frame.positions().begin(), frame.positions().end(), std::back_inserter(oldpos));*/

        //std::ostringstream filename;
        //filename << "mesh_" << step << ".pdb";
        //writePDB(filename.str().c_str(), &pdbOut);

        std::fill(hist.begin(), hist.end(), 0.0);
    }

    int moleculeTypeId = 0;
    for (auto moleculeType : moleculeTypes)
    {
        dr2s[moleculeTypeId] /= moleculeType.count;
        dr2s[moleculeTypeId] /= mds.nsteps();

        std::cout << moleculeType.name << ": " << dr2s[moleculeTypeId] << std::endl;

        moleculeTypeId++;
    }

    /*std::ofstream out;
    out.open("liquid.txt");

    for (int i = 0; i < dr2_liquid.size(); i++)
    {
        double dr2 = dr2_liquid[i];
        double ndr2 = ndr2_liquid[i];
        out << frame.topology()[i].name() << ": " << dr2/ndr2 << " (" << dr2 << ", " << ndr2 << ")" << std::endl;
    }
    out.close();

    out.open("gas.txt");
    for (int i = 0; i < dr2_gas.size(); i++)
    {
        double dr2 = dr2_gas[i];
        double ndr2 = ndr2_gas[i];
        out << frame.topology()[i].name() << ": " << dr2/ndr2 << " (" << dr2 << ", " << ndr2 << ")" << std::endl;
    }
    out.close();*/

/*
    std::cout << std::endl << "Liquid:" << std::endl;
    for (auto& [dr2, ndr2] : zip(dr2_liquid, ndr2_liquid))
    {
        std::cout << dr2/ndr2 << " (" << dr2 << ", " << ndr2 << ")" << std::endl;
    }

    std::cout << std::endl << "Gas:" << std::endl;
    for (auto& [dr2, ndr2] : zip(dr2_gas, ndr2_gas))
    {
        std::cout << dr2/ndr2 << " (" << dr2 << ", " << ndr2 << ")" << std::endl;
    }
*/
    return 0;
}