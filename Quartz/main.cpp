#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>

#include "xyzio.h"

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

#define DIM 3
#define XX 0
#define YY 1
#define ZZ 2

typedef float matrix[DIM][DIM];

typedef float rvec[DIM];

void createBoxMatrix(matrix box, const float a, const float b, const float c, const float alpha, const float beta, const float gamma)
{
    float cosalphastar = (cosf(beta)*cosf(gamma) - cosf(alpha)) / (sinf(beta)*sin(gamma));
    float sinalphastar = sqrtf(1 - cosalphastar*cosalphastar);

    box[XX][XX] = a;
    box[XX][YY] = b*cosf(gamma);
    box[XX][ZZ] = c*cosf(beta);

    box[YY][XX] = 0.0;
    box[YY][YY] = b*sin(gamma);
    box[YY][ZZ] = -c*sin(beta)*cosalphastar;

    box[ZZ][XX] = 0.0;
    box[ZZ][YY] = 0.0;
    box[ZZ][ZZ] = c*sin(beta)*sinalphastar;
}

void calcBoxCenter(const matrix box, rvec center)
{
    center[XX] = 0.0;
    center[YY] = 0.0;
    center[ZZ] = 0.0;
    for (int m = 0; m < DIM; m++)
    {
        for (int d = 0; d < DIM; d++)
        {
            center[d] += 0.5 * box[m][d];
        }
    }
}


void putAtomInTriclinicUC(const matrix box, rvec r, const rvec center)
{
    /* The product of matrix shm with a coordinate gives the shift vector
       which is required determine the periodic cell position */
    float shm01 = box[1][0] / box[1][1];
    float shm02 = (box[1][1] * box[2][0] - box[2][1] * box[1][0]) / (box[1][1] * box[2][2]);
    float shm12 = box[2][1] / box[2][2];

    rvec shiftCenter;
    shiftCenter[XX] = 0.0;
    shiftCenter[YY] = 0.0;
    shiftCenter[ZZ] = 0.0;
    for (int d1 = 0; d1 < DIM; d1++)
    {
        shiftCenter[d1] = 0;
        for (int d2 = 0; d2 < DIM; d2++)
        {
            shiftCenter[d1] += box[d2][d1];
        }
    }
    for (int d = 0; d < DIM; d++)
    {
        shiftCenter[d] *= 0.5;
        shiftCenter[d] -= center[d];
    }

    shiftCenter[0] = shm01 * shiftCenter[1] + shm02 * shiftCenter[2];
    shiftCenter[1] = shm12 * shiftCenter[2];
    shiftCenter[2] = 0;

    for (int m = DIM - 1; m >= 0; m--)
    {
        float shift = shiftCenter[m];
        if (m == 0)
        {
            shift += shm01 * r[1] + shm02 * r[2];
        }
        else if (m == 1)
        {
            shift += shm12 * r[2];
        }
        while (r[m] - shift < 0)
        {
            for (int d = 0; d <= m; d++)
            {
                r[d] += box[m][d];
            }
        }
        while (r[m] - shift >= box[m][m])
        {
            for (int d = 0; d <= m; d++)
            {
                r[d] -= box[m][d];
            }
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc < 8)
    {
        std::cout << "Usags: " << argv[0] << "<input xyz> <input crystal dat> <output file name> <use pbc (yes/no)> <Nx> <Ny> <Lz (in A)>" << std::endl;
    }

    XYZ xyzIn;
    readXYZ(argv[1], &xyzIn);
    
    std::ifstream ucFile;
    ucFile.open(argv[2]);

    std::string outputFilename(argv[3]);

    bool dopbc = false;
    if (strncmp(argv[4], "yes", 3) == 0)
    {
        dopbc = true;
    }

    float a, b, c, alpha, beta, gamma;

    ucFile >> a;
    ucFile >> b;
    ucFile >> c;
    ucFile >> alpha;
    ucFile >> beta;
    ucFile >> gamma;

    std::cout << a << "  " << b << "  " << c << "  " << alpha << "  " << beta << "  " << gamma << std::endl;

    alpha *= M_PI/180.0;
    beta *= M_PI/180.0;
    gamma *= M_PI/180.0;

    int Nx = atoi(argv[5]);
    int Ny = atoi(argv[6]);
    float Lz = atof(argv[7]);
    int Nz = ceil(Lz/c);

    matrix box;
    createBoxMatrix(box, a*Nx, b*Ny, c*Nz, alpha, beta, gamma);

    matrix uc;
    createBoxMatrix(uc, a, b, c, alpha, beta, gamma);

    std::cout << "Recommended x - y box size is " << a*Nx << " - " << b*Ny << std::endl;

    ucFile.close();

    XYZ xyzOut;
    xyzOut.atomCount = xyzIn.atomCount*Nx*Ny;
    xyzOut.atoms = (XYZAtom*)calloc(xyzOut.atomCount, sizeof(XYZAtom));
    for (int iy = 0; iy < Ny; iy ++)
    {
        for (int ix = 0; ix < Nx; ix ++)
        {
            for (int i = 0; i < xyzIn.atomCount; i++)
            {
                int iOut = i + ix*xyzIn.atomCount + iy*Nx*xyzIn.atomCount;
                xyzOut.atoms[iOut].name = xyzIn.atoms[i].name;
                float x = xyzIn.atoms[i].x;
                float y = xyzIn.atoms[i].y;
                float z = xyzIn.atoms[i].z;

                int iz = 0;
                xyzOut.atoms[iOut].x = x + uc[XX][XX]*ix + uc[XX][YY]*iy + uc[XX][ZZ]*iz;
                xyzOut.atoms[iOut].y = y + uc[YY][XX]*ix + uc[YY][YY]*iy + uc[YY][ZZ]*iz;
                xyzOut.atoms[iOut].z = z + uc[ZZ][XX]*ix + uc[ZZ][YY]*iy + uc[ZZ][ZZ]*iz;

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

    for (int iy = 0; iy < Ny; iy ++)
    {
        for (int ix = 0; ix < Nx; ix ++)
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

    FILE* groOut = fopen(outputFilename.c_str(), "w");
    fprintf(groOut, "Quartz slab (%dx%d)\n", Nx, Ny);
    fprintf(groOut, "%d\n", static_cast<int>(atomsOut.size()));
    idx = 0;

    rvec center;
    calcBoxCenter(box, center);
    for (auto atomOut : atomsOut)
    {
        rvec r;
        r[XX] = atomOut.x;
        r[YY] = atomOut.y;
        r[ZZ] = atomOut.z;
        if (dopbc)
        {
            r[XX] -= 0.5*a;
            r[YY] -= 0.5*b;
            r[ZZ] += 0.5*Lz;
            putAtomInTriclinicUC(box, r, center);
        }
        fprintf(groOut, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
            atomOut.ix + Nx*atomOut.iy + 1,
            "Q011",
            atomOut.name.c_str(),
            idx + 1,
            r[XX]*0.1, // A to nm
            r[YY]*0.1,
            r[ZZ]*0.1);
        idx++;
    }
    fprintf(groOut, "%f %f %f %f %f %f %f %f %f\n", 
            box[XX][XX]*Nx*0.1, box[YY][YY]*Ny*0.1, box[ZZ][ZZ]*Nz*0.1,
            box[XX][YY]*Ny*0.1, box[XX][ZZ]*Nz*0.1,
            box[YY][XX]*Nx*0.1, box[YY][ZZ]*Nz*0.1, 
            box[ZZ][XX]*Nx*0.1, box[ZZ][YY]*Ny*0.1);
    fclose(groOut);
}