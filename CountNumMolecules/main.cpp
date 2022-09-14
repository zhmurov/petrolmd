#include <algorithm>
#include <cmath>
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
    float concentration;
    int count;
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
        std::cout << "Unable to open file '" << filename << "'" << std::endl;
    }
    
    std::cout << "Masses database read:" << std::endl;
    for (auto const& atom : masses)
    {
        std::cout << atom.first << ": " << atom.second << std::endl;
    }
}

float getMass(std::string molecule)
{
    std::cout << "Parsing " << molecule << std::endl;
    float totalMass = 0.0f;

    int pos = 0;
    char element;
    int number = 1;

    if (molecule[0] != '#')
    {

        while (pos < molecule.length() && molecule[pos] != ' ' && molecule[pos] != '_')
        {
            char element = molecule[pos];

            std::cout << "Element: " << element << std::endl;
            
            pos++;
            if (molecule[pos] >= '0' && molecule[pos] <= '9')
            {
                int posend = pos + 1;
                while (posend < molecule.length() && (molecule[posend] >= '0' && molecule[posend] <= '9'))
                {
                    posend ++;
                }
                number = atoi(molecule.substr(pos, posend).c_str());
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
    if (argc < 4)
    {
        std::cout << "Usage: ./count_mols <atomic_weights_db.dat> <oil_composition.dat> <output_filename>" << std::endl;
        exit(0);
    }

    readMassesDB(argv[1]);

    char* filename = argv[2];
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
                float massConcentration = atof(tokens[1].c_str())*0.01;

                molecules[tokens[0]].mass = mass;
                molecules[tokens[0]].massConcentration = massConcentration;                

                std::cout << "Mass " << mass << std::endl;
                std::cout << "\% of mass: " << massConcentration << std::endl;
            }
        }
        file.close();

        double totalMassConcentration = 0.0;
        double totalMass = 0.0;
        for (auto const& molecule : molecules)
        {
            std::cout << molecule.first << ": " << molecule.second.mass << ", " << molecule.second.massConcentration << std::endl;
            totalMassConcentration += molecule.second.massConcentration;
            totalMass += molecule.second.mass * molecule.second.massConcentration;
        }
        std::cout << "Total mass concentration: " << totalMassConcentration << std::endl;
        std::cout << "Total mass per 100%: " << totalMass << std::endl;

        double totalConcentration = 0.0;
        for (auto& molecule : molecules)
        {
            molecule.second.concentration = molecule.second.massConcentration / molecule.second.mass;

            molecule.second.count = std::max(static_cast<int>(std::rint(molecule.second.concentration*1000000)), 1);

            std::cout << molecule.first << ": " << molecule.second.concentration  << std::endl;

            totalConcentration += molecule.second.concentration;
        }
        std::cout << "Total concentration: " << totalConcentration << std::endl;

    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
    }

    std::string outputName(argv[3]);
    std::ofstream out;
    out.open(outputName + "_packmol.inp");

    out << "tolerance 2.0" << std::endl;
    out << "filetype pdb" << std::endl;
    out << "output yamburg_1000nm3.pdb" << std::endl;
    out << std::endl;

    for (auto& molecule : molecules)
    {
        out << "structure " << molecule.first << ".pdb" << std::endl;
        out << "  number " << molecule.second.count << std::endl;
        out << "  inside box 0. 0. 0. 1000. 1000. 1000." << std::endl;
        out << "end structure" << std::endl;
    }

    out << std::endl;
    out.close();

    out.open(outputName + ".top");

    out << "; Include forcefield parameters" << std::endl;
    out << "#include \"trappeua.ff/forcefield.itp\"" << std::endl;
    out << "#include \"trappeua.ff/tip4p2005.itp\"" << std::endl;
    for (auto& molecule : molecules)
    {
        out << "#include \"" << molecule.first << ".itp\"" << std::endl;
    }
    out << std::endl;

    out << "[ system ]" << std::endl;
    out << "; Name" << std::endl;
    out << "Yamburg recomb" << std::endl;
    out << std::endl;

    out << "[ molecules ]" << std::endl;
    out << "; Compound     #mols" << std::endl;
    for (auto& molecule : molecules)
    {
        out << molecule.first << "      " << molecule.second.count << std::endl;
    }
    out << std::endl;    

    out.close();
    
    return 0;
}