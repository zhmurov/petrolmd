#include <fstream>
#include <iostream>
#include <string>
#include <iostream>

#include "xyzio.h"

#define NX 2
#define NY 2 

struct Atom {
    std::string name;
    std::string type;
    float charge;
    float x, y, z;
};

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
    xyzOut.atomCount = xyzIn.atomCount*NX*NY;
    xyzOut.atoms = (XYZAtom*)calloc(xyzOut.atomCount, sizeof(XYZAtom));
    for (int ix = 0; ix < NX; ix ++)
    {
        for (int iy = 0; iy < NY; iy ++)
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
    writeXYZ("Quartz/files/slab.xyz", &xyzOut);

    std::ofstream rtpFile;
    rtpFile.open("Quartz/files/sioslab.rtp");

    rtpFile << "[ bondedtypes ]" << std::endl;
    rtpFile << "; Col 1: Type of bond " << std::endl;
    rtpFile << "; Col 2: Type of angles " << std::endl;
    rtpFile << "; Col 3: Type of proper dihedrals " << std::endl;
    rtpFile << "; Col 4: Type of improper dihedrals " << std::endl;
    rtpFile << "; Col 5: Generate all dihedrals if 1, only heavy atoms of 0. " << std::endl;
    rtpFile << "; Col 6: Number of excluded neighbors for nonbonded interactions " << std::endl;
    rtpFile << "; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1 " << std::endl;
    rtpFile << "; Col 8: Remove propers over the same bond as an improper if it is 1 " << std::endl;
    rtpFile << "; bonds  angles  dihedrals  impropers  all_dihedrals  nrexcl  HH14  RemoveDih " << std::endl;
    rtpFile << "1       5        9          2            1           3      1       0" << std::endl;

    rtpFile << "[ SiO_" << NX << "_" << NY << " ]" << std::endl;
    rtpFile << "; " << std::endl;
    rtpFile << "[ atoms ]" << std::endl;
    for (int ix = 0; ix < NX; ix ++)
    {
        for (int iy = 0; iy < NY; iy ++)
        {
            rtpFile << "  SI" << ix << iy << "1     SIH3  0.8000   0" << std::endl;
            rtpFile << "  SI" << ix << iy << "2       SI  1.0800   0" << std::endl;
            rtpFile << "  SI" << ix << iy << "3       SI  1.0800   0" << std::endl;
            rtpFile << "  SI" << ix << iy << "4       SI  1.0800   0" << std::endl;
            rtpFile << "   O" << ix << iy << "5     OSIE -0.5300   0" << std::endl;
            rtpFile << "   O" << ix << iy << "6     OSIE -0.5300   0" << std::endl;
            rtpFile << "   O" << ix << iy << "7     OSIE -0.5300   0" << std::endl;
            rtpFile << "   O" << ix << iy << "8     OSIE -0.5300   0" << std::endl;
            rtpFile << "   O" << ix << iy << "9     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "10     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "11     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "12     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "13     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "14     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "15     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "16     OSIH -0.5400   0" << std::endl;
            rtpFile << "  O" << ix << iy << "17     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "18     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "19     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "20     OSIE -0.5300   0" << std::endl;
            rtpFile << " SI" << ix << iy << "21       SI  1.0800   0" << std::endl;
            rtpFile << " SI" << ix << iy << "22       SI  1.0800   0" << std::endl;
            rtpFile << " SI" << ix << iy << "23       SI  1.0800   0" << std::endl;
            rtpFile << " SI" << ix << iy << "24       SI  1.0800   0" << std::endl;
            rtpFile << " SI" << ix << iy << "25     SIH3  0.8000   0" << std::endl;
            rtpFile << " SI" << ix << iy << "26       SI  1.0800   0" << std::endl;
            rtpFile << " SI" << ix << iy << "27       SI  1.0800   0" << std::endl;
            rtpFile << " SI" << ix << iy << "28       SI  1.0800   0" << std::endl;
            rtpFile << "  O" << ix << iy << "29     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "30     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "31     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "32     OSIH -0.5400   0" << std::endl;
            rtpFile << "  O" << ix << iy << "33     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "34     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "35     OSIE -0.5300   0" << std::endl;
            rtpFile << "  O" << ix << iy << "36     OSIE -0.5300   0" << std::endl;
            rtpFile << "  H" << ix << iy << "37     HSIA -0.1500   0" << std::endl;
            rtpFile << "  H" << ix << iy << "38     HSIO  0.3200   0" << std::endl;
            rtpFile << "  H" << ix << iy << "39     HSIA -0.1500   0" << std::endl;
            rtpFile << "  H" << ix << iy << "40     HSIO  0.3200   0" << std::endl;
        }
    }
    rtpFile << "[ bonds ]" << std::endl;
    for (int ix = 0; ix < NX; ix ++)
    {
        for (int iy = 0; iy < NY; iy ++)
        {
            rtpFile << "   O" << ix << iy << "5   SI" << ix << iy << "1" << std::endl;
            rtpFile << "   O" << ix << iy << "9   SI" << ix << iy << "1" << std::endl;
            rtpFile << "  H" << ix << iy << "37   SI" << ix << iy << "1" << std::endl;
            rtpFile << "   O" << ix << iy << "6   SI" << ix << iy << "2" << std::endl;
            rtpFile << "  O" << ix << iy << "10   SI" << ix << iy << "2" << std::endl;
            rtpFile << "  O" << ix << iy << "13   SI" << ix << iy << "2" << std::endl;
            rtpFile << "  O" << ix << iy << "18   SI" << ix << iy << "2" << std::endl;
            rtpFile << "   O" << ix << iy << "7   SI" << ix << iy << "3" << std::endl;
            rtpFile << "  O" << ix << iy << "11   SI" << ix << iy << "3" << std::endl;
            rtpFile << "  O" << ix << iy << "14   SI" << ix << iy << "3" << std::endl;
            rtpFile << "  O" << ix << iy << "19   SI" << ix << iy << "3" << std::endl;
            rtpFile << "  O" << ix << iy << "15   SI" << ix << iy << "4" << std::endl;
            rtpFile << "  O" << ix << iy << "20   SI" << ix << iy << "4" << std::endl;
            rtpFile << " SI" << ix << iy << "21    O" << ix << iy << "5" << std::endl;
            rtpFile << " SI" << ix << iy << "24    O" << ix << iy << "8" << std::endl;
            rtpFile << " SI" << ix << iy << "26   O" << ix << iy << "10" << std::endl;
            rtpFile << " SI" << ix << iy << "27   O" << ix << iy << "11" << std::endl;
            rtpFile << " SI" << ix << iy << "28   O" << ix << iy << "12" << std::endl;
            rtpFile << " SI" << ix << iy << "25   O" << ix << iy << "13" << std::endl;
            rtpFile << " SI" << ix << iy << "26   O" << ix << iy << "14" << std::endl;
            rtpFile << " SI" << ix << iy << "27   O" << ix << iy << "15" << std::endl;
            rtpFile << " SI" << ix << iy << "28   O" << ix << iy << "16" << std::endl;
            rtpFile << " SI" << ix << iy << "21   O" << ix << iy << "17" << std::endl;
            rtpFile << " SI" << ix << iy << "23   O" << ix << iy << "19" << std::endl;
            rtpFile << " SI" << ix << iy << "24   O" << ix << iy << "20" << std::endl;
            rtpFile << "  O" << ix << iy << "29  SI" << ix << iy << "21" << std::endl;
            rtpFile << "  O" << ix << iy << "33  SI" << ix << iy << "21" << std::endl;
            rtpFile << "  O" << ix << iy << "30  SI" << ix << iy << "22" << std::endl;
            rtpFile << "  O" << ix << iy << "34  SI" << ix << iy << "22" << std::endl;
            rtpFile << "  O" << ix << iy << "31  SI" << ix << iy << "23" << std::endl;
            rtpFile << "  O" << ix << iy << "35  SI" << ix << iy << "23" << std::endl;
            rtpFile << "  O" << ix << iy << "32  SI" << ix << iy << "24" << std::endl;
            rtpFile << "  O" << ix << iy << "36  SI" << ix << iy << "24" << std::endl;
            rtpFile << "  H" << ix << iy << "39  SI" << ix << iy << "25" << std::endl;
            rtpFile << "  O" << ix << iy << "29  SI" << ix << iy << "26" << std::endl;
            rtpFile << "  O" << ix << iy << "34  SI" << ix << iy << "26" << std::endl;
            rtpFile << "  O" << ix << iy << "36  SI" << ix << iy << "28" << std::endl;
            rtpFile << "  H" << ix << iy << "40   O" << ix << iy << "32" << std::endl;
            rtpFile << "  O" << ix << iy << "16   H" << ix << iy << "38" << std::endl;
        }
    }
    rtpFile.close();

}