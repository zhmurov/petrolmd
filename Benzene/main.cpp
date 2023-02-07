#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#define CC_BOND_DISTANCE 0.1375
#define CH_BOND_DISTANCE 0.1080

void addAtom(FILE* out, const char* atomName, const int id, const double x, const double y, const double z)
{
    fprintf(out, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "BENZ", atomName, id, x, y, z);
}

int main(int argc, char *argv[])
{
    double a = CC_BOND_DISTANCE;
    double RC = a;
    double rC = a*cos(M_PI/6.0);
    double apb = CC_BOND_DISTANCE + CH_BOND_DISTANCE;
    double RH = apb;
    double rH = apb*cos(M_PI/6.0);


    FILE* groOut = fopen("C6H6.gro", "w");
    fprintf(groOut, "Benzene\n");
    fprintf(groOut, "12\n");
    addAtom(groOut, "CG", 1, RC, 0.0, 0.0);
    addAtom(groOut, "HG", 2, RH, 0.0, 0.0);
    addAtom(groOut, "CD1", 3, 0.5*RC, rC, 0.0);
    addAtom(groOut, "HD1", 4, 0.5*RH, rH, 0.0);
    addAtom(groOut, "CD2", 5, 0.5*RC, -rC, 0.0);
    addAtom(groOut, "HD2", 6, 0.5*RH, -rH, 0.0);
    addAtom(groOut, "CE1", 7, -0.5*RC, rC, 0.0);
    addAtom(groOut, "HE1", 8, -0.5*RH, rH, 0.0);
    addAtom(groOut, "CE2", 9, -0.5*RC, -rC, 0.0);
    addAtom(groOut, "HE2", 10, -0.5*RH, -rH, 0.0);
    addAtom(groOut, "CZ", 11, -RC, 0.0, 0.0);
    addAtom(groOut, "HZ", 12, -RH, 0.0, 0.0);
    fclose(groOut);
}
