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

struct Coord
{
    float x;
    float y;
    float z;
};

struct Atom
{
    std::string name;
    std::string type;
    float charge;
    int cgr;
    Coord r;
    int ix;
    int iy;
};

#define DIM 3
#define XX 0
#define YY 1
#define ZZ 2

typedef float matrix[DIM][DIM];

double cot(double alpha)
{
    return 1.0/tan(alpha);
}

double csc(double alpha)
{
    return 1.0/sin(alpha);
}

void createBoxMatrix(matrix box, const float a, const float b, const float c, const float alpha, const float beta, const float gamma)
{
    float cosalphastar = (cos(beta)*cos(gamma) - cos(alpha)) / (sin(beta)*sin(gamma));
    float sinalphastar = sqrt(1 - cosalphastar*cosalphastar);

    box[XX][XX] = a;
    box[XX][YY] = b*cos(gamma);
    box[XX][ZZ] = c*cos(beta);

    box[YY][XX] = 0.0;
    box[YY][YY] = b*sin(gamma);
    box[YY][ZZ] = -c*sin(beta)*cosalphastar;

    box[ZZ][XX] = 0.0;
    box[ZZ][YY] = 0.0;
    box[ZZ][ZZ] = c*sin(beta)*sinalphastar;
}

void createInverseBoxMatrix(matrix ibox, const float a, const float b, const float c, const float alpha, const float beta, const float gamma)
{
    float cosalphastar = (cosf(beta)*cosf(gamma) - cosf(alpha)) / (sinf(beta)*sin(gamma));
    float sinalphastar = sqrtf(1 - cosalphastar*cosalphastar);
    float cscalphastar = 1.0/sinalphastar;
    float cotalphastar = cosalphastar/sinalphastar;

    ibox[XX][XX] = 1.0/a;
    ibox[XX][YY] = -cot(gamma)/a;
    ibox[XX][ZZ] = -(cscalphastar*(cot(beta) + cot(gamma)*cosalphastar))/a;

    ibox[YY][XX] = 0.0;
    ibox[YY][YY] = csc(gamma)/b;
    ibox[YY][ZZ] = csc(gamma)*cotalphastar/b;

    ibox[ZZ][XX] = 0.0;
    ibox[ZZ][YY] = 0.0;
    ibox[ZZ][ZZ] = csc(beta)*cscalphastar/c;

}

Coord calcBoxCenter(const matrix box)
{
    Coord center;
    center.x = 0.0;
    center.y = 0.0;
    center.z = 0.0;
    for (int d = 0; d < DIM; d++)
    {
        center.x += 0.5 * box[d][XX];
        center.y += 0.5 * box[d][YY];
        center.z += 0.5 * box[d][ZZ];
    }
    return center;
}

Coord putAtomInTriclinicBox(const matrix box, Coord position, const Coord center)
{
    /* The product of matrix shm with a coordinate gives the shift vector
       which is required determine the periodic cell position */
    float shm01 = box[1][0] / box[1][1];
    float shm02 = (box[1][1] * box[2][0] - box[2][1] * box[1][0]) / (box[1][1] * box[2][2]);
    float shm12 = box[2][1] / box[2][2];

    Coord shiftCenter;
    shiftCenter.x = 0.0;
    shiftCenter.y = 0.0;
    shiftCenter.z = 0.0;
    for (int d = 0; d < DIM; d++)
    {
        shiftCenter.x += box[d][XX];
        shiftCenter.y += box[d][YY];
        shiftCenter.z += box[d][ZZ];
    }
    shiftCenter.x *= 0.5;
    shiftCenter.y *= 0.5;
    shiftCenter.z *= 0.5;
    shiftCenter.x -= center.x;
    shiftCenter.y -= center.y;
    shiftCenter.z -= center.z;
    
    shiftCenter.x = shm01 * shiftCenter.y + shm02 * shiftCenter.z;
    shiftCenter.y = shm12 * shiftCenter.z;
    shiftCenter.z = 0;

    Coord r = position;

    while (r.z - shiftCenter.z < 0)
    {
        r.x += box[ZZ][XX];
        r.y += box[ZZ][YY];
        r.z += box[ZZ][ZZ];
    }
    while (r.z - shiftCenter.z >= box[ZZ][ZZ])
    {
        r.x -= box[ZZ][XX];
        r.y -= box[ZZ][YY];
        r.z -= box[ZZ][ZZ];
    }
    shiftCenter.y += shm12 * r.z;
     
    while (r.y - shiftCenter.y < 0)
    {
        r.x += box[YY][XX];
        r.y += box[YY][YY];
        r.z += box[YY][ZZ];
    }
    while (r.y - shiftCenter.y >= box[YY][YY])
    {
        r.x -= box[YY][XX];
        r.y -= box[YY][YY];
        r.z -= box[YY][ZZ];
    }

    shiftCenter.x += shm01 * r.y + shm02 * r.z;
    while (r.x - shiftCenter.x < 0)
    {
        r.x += box[XX][XX];
        r.y += box[XX][YY];
        r.z += box[XX][ZZ];
    }
    while (r.x - shiftCenter.x >= box[XX][XX])
    {
        r.x -= box[XX][XX];
        r.y -= box[XX][YY];
        r.z -= box[XX][ZZ];
    }

    return r;
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
        atomsOut[i].r.x = xyzOut.atoms[i].x;
        atomsOut[i].r.y = xyzOut.atoms[i].y;
        atomsOut[i].r.z = xyzOut.atoms[i].z;
    }

    FILE* groOut = fopen(outputFilename.c_str(), "w");
    fprintf(groOut, "Quartz slab (%dx%d)\n", Nx, Ny);
    fprintf(groOut, "%d\n", static_cast<int>(atomsOut.size()));
    idx = 0;

    Coord center = calcBoxCenter(box);
    for (auto atomOut : atomsOut)
    {
        Coord r = atomOut.r;
        if (dopbc)
        {
            /*r.x -= 0.5*a;
            r.y -= 0.5*b;
            r.z += 0.5*Lz;*/
            r = putAtomInTriclinicBox(box, r, center);
        }
        fprintf(groOut, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
            atomOut.ix + Nx*atomOut.iy + 1,
            "Q011",
            atomOut.name.c_str(),
            idx + 1,
            r.x*0.1, // A to nm
            r.y*0.1,
            r.z*0.1);
        idx++;
    }
    fprintf(groOut, "%f %f %f %f %f %f %f %f %f\n", 
            box[XX][XX]*Nx*0.1, box[YY][YY]*Ny*0.1, box[ZZ][ZZ]*Nz*0.1,
            box[XX][YY]*Ny*0.1, box[XX][ZZ]*Nz*0.1,
            box[YY][XX]*Nx*0.1, box[YY][ZZ]*Nz*0.1, 
            box[ZZ][XX]*Nx*0.1, box[ZZ][YY]*Ny*0.1);
    fclose(groOut);
}