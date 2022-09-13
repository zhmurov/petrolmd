#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <sstream>
#include <vector>

std::map<char, float> masses;

struct Molecule
{
    float mass;
    float massConcentration;
};

std::map<std::string, Molecule> molecules;

#define Lx 100.0f
#define Ly 100.0f
#define Lz 100.0f

void readMassesDB(std::string filename)
{
    std::string line;
    std::ifstream file(filename);
    if (file.is_open())
    {
        while (getline(file, line))
        {
            char element = line.c_str()[0];

            if (element != '#')
            {
                std::string::size_type pos = line.find_first_of("1234567890");
                std::string mass = line.substr(pos, line.length());
                masses[element] = atof(mass.c_str());
            }
        }
        file.close();
    }
    else
    {
        std::cout << "Unable to open file '" << filename << std::endl;
    }
    
    std::cout << "Masses database read:" << std::endl;
    for (auto const& atom : masses)
    {
        std::cout << atom.first << ": " << atom.second << std::endl;
    }
}

float getMass(std::string compound)
{
    std::cout << "Parsing " << compound << std::endl;
    float totalMass = 0.0f;

    int pos = 0;
    char element;
    int number = 1;

    if (compound[0] != '#')
    {

        while (pos < compound.length() && compound[pos] != ' ' && compound[pos] != '_')
        {
            char element = compound[pos];

            std::cout << "Element: " << element << std::endl;
            
            pos++;
            if (compound[pos] >= '0' && compound[pos] <= '9')
            {
                int posend = pos + 1;
                while (posend < compound.length() && (compound[posend] >= '0' && compound[posend] <= '9'))
                {
                    posend ++;
                }
                number = atoi(compound.substr(pos, posend).c_str());
                pos = posend;
            }
            else
            {
                number = 1;
            }
            totalMass += masses[element]*number;
        }
    }

    return totalMass;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Please, specify the data file with the list of compounds." << std::endl;
        exit(0);
    }

    readMassesDB("atomic_weights.dat");

    char* filename = argv[1];
    std::cout << "Reading: " << filename << "\n";
    
    std::string line;
    std::ifstream file(filename);
    if (file.is_open())
    {
        while (getline(file, line))
        {
            std::istringstream iss(line);
            std::vector<std::string> tokens(std::istream_iterator<std::string>{iss},
                                 std::istream_iterator<std::string>());

            //std::cout << "Mass: " << tokens[0] << "  Concentration: " << tokens[1] << std::endl;

            if (tokens[0][0] != '#')
            {
                float mass = getMass(tokens[0]);
                float massConcentration = atof(tokens[1].c_str());

                molecules[tokens[0]].mass = mass;
                molecules[tokens[0]].massConcentration = massConcentration;                

                std::cout << "Mass " << mass << std::endl;
                std::cout << "\% of mass: " << massConcentration << std::endl;
            }
        }
        file.close();
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
    }
    
    return 0;
}