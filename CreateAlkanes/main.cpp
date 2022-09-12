#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <string>

#include "pdbio.h"

#define CC_BOND_DISTANCE 1.540
#define CCC_BOND_ANGLE   114.0
#define MAX_ATOM_COUNT   100

int main(int argc, char *argv[])
{
    PDB alkanePDB;
    alkanePDB.atoms = (PDBAtom*)calloc(MAX_ATOM_COUNT, sizeof(PDBAtom));

    for (int atomCount = 1; atomCount <= MAX_ATOM_COUNT; atomCount++)
    {

        alkanePDB.atomCount = atomCount;
        
        float x = 0.0;
        float y = 0.0;
        float z = 0.0;

        float d = CC_BOND_DISTANCE;
        float alpha = (0.5*CCC_BOND_ANGLE)*M_PI/180.0;

        for (int i = 0; i < atomCount; i++)
        {
            alkanePDB.atoms[i].id = i + 1;
            sprintf(alkanePDB.atoms[i].name, "C");
            alkanePDB.atoms[i].chain = ' ';
            std::string resName = (i == 0 || i == atomCount - 1) ? "CH3" : "CH2";
            sprintf(alkanePDB.atoms[i].resName, "%s", resName.c_str());

            alkanePDB.atoms[i].altLoc = ' ';
            alkanePDB.atoms[i].resid = i + 1;
            alkanePDB.atoms[i].x = x;
            alkanePDB.atoms[i].y = y;
            alkanePDB.atoms[i].z = z;

            sprintf(alkanePDB.atoms[i].segment, " ");

            alkanePDB.atoms[i].occupancy = 0.0;
            alkanePDB.atoms[i].beta = 0.0;

            x += d*sinf(alpha);
            y = (i % 2) ? 0.0 : d*cosf(alpha);
        }

        std::string filename = "C" + std::to_string(atomCount) + "H" + std::to_string(atomCount*2 + 2) + ".pdb";

        writePDB(filename.c_str(), &alkanePDB);
    }

}