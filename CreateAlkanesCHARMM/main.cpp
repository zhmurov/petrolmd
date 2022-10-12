#include <string>
#include <iostream>
#include <fstream>

#define MAX_ATOM_COUNT   20

int main(int argc, char *argv[])
{
    for (int atomCount = 1; atomCount <= MAX_ATOM_COUNT; atomCount++)
    {
        std::ofstream file;
        std::string molname = "C" + std::to_string(atomCount) + "H" + std::to_string(atomCount*2 + 2);
        file.open(molname + ".top");

        file << "[ " << molname << " ]" << std::endl;
        file << "  [ atoms ]" << std::endl;
        
        file << "       C1    CC33A -0.2700   1" << std::endl;
        file << "      H11    HCA3A  0.0900   1" << std::endl;
        file << "      H12    HCA3A  0.0900   1" << std::endl;
        file << "      H13    HCA3A  0.0900   1" << std::endl;

        for (int i = 1; i < atomCount - 1; i++)
        {
            file << "       C" << i + 1 << "    CC32A -0.1800   " << i + 1 << std::endl;
            file << "      H" << i + 1 << "1    HCA2A  0.0900   " << i + 1 << std::endl;
            file << "      H" << i + 1 << "2    HCA2A  0.0900   " << i + 1 << std::endl;
        }

        file << "       C" << atomCount << "    CC33A -0.2700   " << atomCount << std::endl;
        file << "      H" << atomCount << "1    HCA3A  0.0900   " << atomCount << std::endl;
        file << "      H" << atomCount << "2    HCA3A  0.0900   " << atomCount << std::endl;
        file << "      H" << atomCount << "3    HCA3A  0.0900   " << atomCount << std::endl;

        file << "  [ bonds ]" << std::endl;
        for (int i = 1; i < atomCount; i++)
        {
            file << "      C" << i << "    C" << i + 1 << std::endl; 
        }

        file << "      C1    H11" << std::endl;
        file << "      C1    H12" << std::endl;
        file << "      C1    H13" << std::endl;
       
        for (int i = 2; i < atomCount; i++)
        {
            file << "      C" << i << "    H" << i << "1" << std::endl;
            file << "      C" << i << "    H" << i << "2" << std::endl;
        }

        file << "      C" << atomCount << "    H" << atomCount << "1" << std::endl;
        file << "      C" << atomCount << "    H" << atomCount << "2" << std::endl;
        file << "      C" << atomCount << "    H" << atomCount << "3" << std::endl;

        file.close();
    }

}