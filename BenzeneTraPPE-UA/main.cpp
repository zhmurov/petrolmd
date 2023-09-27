#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#define CC_BOND_DISTANCE 0.14

void addAtom(FILE* out, const char* atomName, const int id, const double x, const double y, const double z)
{
    fprintf(out, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "BENZ", atomName, id, x, y, z);
}

int main(int argc, char *argv[])
{
    double a = CC_BOND_DISTANCE;
    double RC = a;
    double rC = a*cos(M_PI/6.0);

    FILE* groOut = fopen("C6H6.gro", "w");
    fprintf(groOut, "Benzene\n");
    fprintf(groOut, "6\n");
    addAtom(groOut, "C1", 1, -RC, 0.0, 0.0);
    addAtom(groOut, "C2", 2, -0.5*RC, rC, 0.0);
    addAtom(groOut, "C3", 3, 0.5*RC, rC, 0.0);    
    addAtom(groOut, "C4", 4, RC, 0.0, 0.0);
    addAtom(groOut, "C5", 5, 0.5*RC, -rC, 0.0);
    addAtom(groOut, "C6", 6, -0.5*RC, -rC, 0.0);
  
    fclose(groOut);
}
