#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "pdbio.h"

#define CC_BOND_DISTANCE 1.530
#define CCC_BOND_ANGLE   114.0
#define CH_BOND_DISTANCE 1.111
#define HCH_BOND_ANGLE   109.0
#define MAX_ATOM_COUNT   100

#define CH_BOND_DISTANCE_METHANE 1.087

int main(int argc, char *argv[])
{
    PDB alkanePDB;
    alkanePDB.atoms = (PDBAtom*)calloc(MAX_ATOM_COUNT*4, sizeof(PDBAtom));

    for (int atomCount = 1; atomCount <= MAX_ATOM_COUNT; atomCount++)
    {    
        std::string molName = "C" + std::to_string(atomCount) + "H" + std::to_string(atomCount*2 + 2);

        float x = 0.0;
        float y = 0.0;
        float z = 0.0;

        float d = CC_BOND_DISTANCE;
        float alpha = (0.5*CCC_BOND_ANGLE)*M_PI/180.0;

        int atomIndex = 0;

        for (int i = 0; i < atomCount; i++)
        {
            alkanePDB.atoms[atomIndex].id = i + 1;
            sprintf(alkanePDB.atoms[atomIndex].name, "C");
            alkanePDB.atoms[atomIndex].chain = ' ';

            alkanePDB.atoms[atomIndex].altLoc = ' ';
            alkanePDB.atoms[atomIndex].resid = i + 1;
            alkanePDB.atoms[atomIndex].x = x;
            alkanePDB.atoms[atomIndex].y = y;
            alkanePDB.atoms[atomIndex].z = z;

            sprintf(alkanePDB.atoms[atomIndex].segment, " ");

            alkanePDB.atoms[atomIndex].occupancy = 0.0;
            alkanePDB.atoms[atomIndex].beta = 0.0;

            int numHydrogens;
            std::string resName;
            if (atomCount == 1)
            {
                numHydrogens = 4;
                resName = "CH4";
            }
            else if (i == 0 || i == atomCount - 1)
            {
                numHydrogens = 3;
                resName = "CH3";
            } 
            else
            {
                numHydrogens = 2;
                resName = "CH2";
            }
            sprintf(alkanePDB.atoms[atomIndex].resName, "%s", resName.c_str());
            for (int deltaIndex = 1; deltaIndex <= numHydrogens; deltaIndex ++)
            {
                alkanePDB.atoms[atomIndex + deltaIndex].id = i + 1;
                alkanePDB.atoms[atomIndex + deltaIndex].chain = ' ';
                sprintf(alkanePDB.atoms[atomIndex + deltaIndex].resName, "%s", resName.c_str());

                alkanePDB.atoms[atomIndex + deltaIndex].altLoc = ' ';
                alkanePDB.atoms[atomIndex + deltaIndex].resid = i + 1;

                sprintf(alkanePDB.atoms[atomIndex + deltaIndex].segment, " ");

                alkanePDB.atoms[atomIndex + deltaIndex].occupancy = 0.0;
                alkanePDB.atoms[atomIndex + deltaIndex].beta = 0.0;
            }

            float dH = CH_BOND_DISTANCE;
            float alphaH = (0.5*HCH_BOND_ANGLE)*M_PI/180.0; 

            float dx, dy, dy2, dz;

            // The methane molecule has 4 hydrogens and dealt with separately (below)
            if (atomCount != 1)
            {
                // First and last carbons have three hydrogens
                if (i == 0 || i == atomCount - 1)
                {
                    dx = dH*sqrtf(1.0-(4.0/3.0)*sinf(alphaH)*sinf(alphaH));
                    dy = dH*2.0*sinf(alphaH)/sqrtf(3);
                    dy2 = dH*sinf(alphaH)/sqrtf(3);
                    dz = dH*sinf(alphaH);
                    if (i == atomCount - 1 && atomCount % 2 == 0)
                    {
                        dy = -dy;
                        dy2 = -dy2;
                    }
                }
                else
                {
                    dx = 0.0f;
                    dy = dH*cosf(alphaH);
                    dz = dH*sinf(alphaH);
                    dy2 = 0.0f;
                }

                // First carbon
                if (i == 0)
                {
                    sprintf(alkanePDB.atoms[atomIndex + 1].name, "H1");
                    alkanePDB.atoms[atomIndex + 1].x = x - dx;
                    alkanePDB.atoms[atomIndex + 1].y = y + dy;
                    alkanePDB.atoms[atomIndex + 1].z = z;

                    sprintf(alkanePDB.atoms[atomIndex + 2].name, "H2");
                    alkanePDB.atoms[atomIndex + 2].x = x - dx;
                    alkanePDB.atoms[atomIndex + 2].y = y - dy2;
                    alkanePDB.atoms[atomIndex + 2].z = z + dz;

                    sprintf(alkanePDB.atoms[atomIndex + 3].name, "H3");
                    alkanePDB.atoms[atomIndex + 3].x = x - dx;
                    alkanePDB.atoms[atomIndex + 3].y = y - dy2;
                    alkanePDB.atoms[atomIndex + 3].z = z - dz;

                    atomIndex += 3;
                }
                // Last carbon
                else if (i == atomCount - 1)
                {
                    sprintf(alkanePDB.atoms[atomIndex + 1].name, "H1");
                    alkanePDB.atoms[atomIndex + 1].x = x + dx;
                    alkanePDB.atoms[atomIndex + 1].y = y + dy;
                    alkanePDB.atoms[atomIndex + 1].z = z;

                    sprintf(alkanePDB.atoms[atomIndex + 2].name, "H2");
                    alkanePDB.atoms[atomIndex + 2].x = x + dx;
                    alkanePDB.atoms[atomIndex + 2].y = y - dy2;
                    alkanePDB.atoms[atomIndex + 2].z = z + dz;

                    sprintf(alkanePDB.atoms[atomIndex + 3].name, "H3");
                    alkanePDB.atoms[atomIndex + 3].x = x + dx;
                    alkanePDB.atoms[atomIndex + 3].y = y - dy2;
                    alkanePDB.atoms[atomIndex + 3].z = z - dz;

                    atomIndex += 3;
                }
                // Carbon atoms in the middle
                else
                {
                    sprintf(alkanePDB.atoms[atomIndex + 1].name, "H1");
                    alkanePDB.atoms[atomIndex + 1].x = x;
                    alkanePDB.atoms[atomIndex + 1].y = (i % 2) ? y + dy : y - dy;
                    alkanePDB.atoms[atomIndex + 1].z = z + dz;

                    sprintf(alkanePDB.atoms[atomIndex + 2].name, "H2");
                    alkanePDB.atoms[atomIndex + 2].x = x;
                    alkanePDB.atoms[atomIndex + 2].y = (i % 2) ? y + dy : y - dy;
                    alkanePDB.atoms[atomIndex + 2].z = z - dz;

                    atomIndex += 2;
                }
            }
            // Methane (4 hydrogens)
            else
            {
                float r = CH_BOND_DISTANCE_METHANE;

                sprintf(alkanePDB.atoms[atomIndex + 1].name, "H1");
                alkanePDB.atoms[atomIndex + 1].x = x + r * sqrtf(8.0/9.0);
                alkanePDB.atoms[atomIndex + 1].y = y + r * 0.0;
                alkanePDB.atoms[atomIndex + 1].z = z - r * 1.0/3.0;

                sprintf(alkanePDB.atoms[atomIndex + 2].name, "H2");
                alkanePDB.atoms[atomIndex + 2].x = x - r * sqrtf(2.0/9.0);
                alkanePDB.atoms[atomIndex + 2].y = y + r * sqrtf(2.0/3.0);
                alkanePDB.atoms[atomIndex + 2].z = z - r * 1.0/3.0;

                sprintf(alkanePDB.atoms[atomIndex + 3].name, "H3");
                alkanePDB.atoms[atomIndex + 3].x = x - r * sqrtf(2.0/9.0);
                alkanePDB.atoms[atomIndex + 3].y = y - r * sqrtf(2.0/3.0);
                alkanePDB.atoms[atomIndex + 3].z = z - r * 1.0/3.0;

                sprintf(alkanePDB.atoms[atomIndex + 4].name, "H4");
                alkanePDB.atoms[atomIndex + 4].x = x + r * 0.0;
                alkanePDB.atoms[atomIndex + 4].y = y + r * 0.0;
                alkanePDB.atoms[atomIndex + 4].z = z + r * 1.0;

                atomIndex += 4;
            }
            
            x += d*sinf(alpha);
            y = (i % 2) ? 0.0 : d*cosf(alpha);
            atomIndex ++;
        }

        std::string filename = molName + ".pdb";

        alkanePDB.atomCount = atomIndex;

        writePDB(filename.c_str(), &alkanePDB);

    }  

}