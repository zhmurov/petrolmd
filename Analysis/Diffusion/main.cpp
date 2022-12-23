#include <iostream>
#include <sstream>
#include <tuple>
#include <ranges>
#include <chemfiles.hpp>

#include "pdbio.h"

#define NBINX 50
#define NBINY 50
#define NBINZ 50

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


int main() {
    int nbin[3];
    nbin[0] = NBINX;
    nbin[1] = NBINY;
    nbin[2] = NBINZ;

    //chemfiles::Trajectory mds("/home/zhmurov/Data/Sirius/PetrolMD/25-charmm-yamburg/md_iso.xtc");
    //mds.set_topology("/home/zhmurov/Data/Sirius/PetrolMD/25-charmm-yamburg/md_iso.tpr");

    chemfiles::Trajectory mds("/home/zhmurov/Sirius/5-charmm/md_iso.xtc");
    mds.set_topology("/home/zhmurov/Sirius/5-charmm/md_iso.tpr");

    std::vector<double> hist = std::vector<double>(NBINX*NBINY*NBINZ, 0.0);

    PDB pdbOut;
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
    }
    auto frame = mds.read();
    
    std::vector<chemfiles::Vector3D> oldpos = std::vector<chemfiles::Vector3D>(frame.size(), {0.0, 0.0, 0.0});
    std::vector<double> dr2_liquid = std::vector<double>(frame.size(), 0.0);
    std::vector<double> dr2_gas = std::vector<double>(frame.size(), 0.0);
    std::vector<int> ndr2_liquid = std::vector<int>(frame.size(), 0);
    std::vector<int> ndr2_gas = std::vector<int>(frame.size(), 0);

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

        for (int ix = 0; ix < NBINX; ix++)
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
        }

        int i = 0;
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
            if (hist[getIdx(idx)] > 3)
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

        std::copy(frame.positions().begin(), frame.positions().end(), std::back_inserter(oldpos));

        //std::ostringstream filename;
        //filename << "mesh_" << step << ".pdb";
        //writePDB(filename.str().c_str(), &pdbOut);

        std::fill(hist.begin(), hist.end(), 0.0);
    }

    std::cout << std::endl << "Liquid:" << std::endl;
    for (int i = 0; i < dr2_liquid.size(); i++)
    {
        double dr2 = dr2_liquid[i];
        double ndr2 = ndr2_liquid[i];
        std::cout << frame.topology()[i].name() << ": " << dr2/ndr2 << " (" << dr2 << ", " << ndr2 << ")" << std::endl;
    }

    std::cout << std::endl << "Gas:" << std::endl;
    for (int i = 0; i < dr2_gas.size(); i++)
    {
        double dr2 = dr2_gas[i];
        double ndr2 = ndr2_gas[i];
        std::cout << frame.topology()[i].name() << ": " << dr2/ndr2 << " (" << dr2 << ", " << ndr2 << ")" << std::endl;
    }

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