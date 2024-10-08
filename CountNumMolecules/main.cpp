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
    std::string name;
    float mass;
    float massConcentration;
    float concentration;
    int count;

    float xmin, ymin, zmin, xmax, ymax, zmax;
};

std::vector<Molecule> molecules;

float Lx = 10.0f;
float Ly = 10.0f;
float Lz = 10.0f;

float density = 1000.0f; // g/l

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
    if (argc < 7)
    {
        std::cout << "Usage: ./count_mols <atomic_weights_db.dat> <oil_composition.dat> <output_filename> <Lx(nm)> <Ly(nm)> <Lz(nm)>" << std::endl;
        exit(0);
    }

    Lx = atof(argv[4]);
    Ly = atof(argv[5]);
    Lz = atof(argv[6]);

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
                float massConcentration = std::stof(tokens[1])*0.01;

                Molecule molecule;
                molecule.name = tokens[0];
                molecule.mass = mass;
                molecule.massConcentration = massConcentration;

                if (tokens.size() == 8)
                {
                    molecule.xmin = std::stof(tokens[2]);
                    molecule.ymin = std::stof(tokens[3]);
                    molecule.zmin = std::stof(tokens[4]);

                    molecule.xmax = std::stof(tokens[5]);
                    molecule.ymax = std::stof(tokens[6]);
                    molecule.zmax = std::stof(tokens[7]);
                }
                else
                {
                    molecule.xmin = 0.0;
                    molecule.ymin = 0.0;
                    molecule.zmin = 0.0;

                    molecule.xmax = Lx;
                    molecule.ymax = Ly;
                    molecule.zmax = Lz;
                }

                molecules.push_back(molecule);             

                std::cout << "Mass " << mass << std::endl;
                std::cout << "\% of mass: " << massConcentration << std::endl;
            }
        }
        file.close();

        double totalMassConcentration = 0.0;
        double totalMass = 0.0;
        for (auto const& molecule : molecules)
        {
            std::cout << molecule.name << ": " << molecule.mass << ", " << molecule.massConcentration << std::endl;
            totalMassConcentration += molecule.massConcentration;
            totalMass += molecule.mass * molecule.massConcentration;
        }
        std::cout << "Total mass concentration: " << totalMassConcentration << std::endl;
        std::cout << "Total mass per 100%: " << totalMass << std::endl;

        double alpha = (density*0.6022*Lx*Ly*Lz) / totalMassConcentration;

        std::cout << "Density scaling factor: " << alpha << std::endl;

        double totalConcentration = 0.0;
        for (auto& molecule : molecules)
        {
            molecule.concentration = alpha * molecule.massConcentration / molecule.mass;

            molecule.count = std::max(static_cast<int>(std::rint(molecule.concentration)), 1);

            std::cout << molecule.name << ": " << molecule.concentration  << std::endl;

            totalConcentration += molecule.concentration;
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
    out << "output " << outputName << ".pdb" << std::endl;
    out << std::endl;

    for (auto& molecule : molecules)
    {
        out << "structure toppar/" << molecule.name << ".pdb" << std::endl;
        out << "  number " << molecule.count << std::endl;
        out << "  inside box " 
            << molecule.xmin*10.0 << " "
            << molecule.ymin*10.0 << " "
            << molecule.zmin*10.0 << " "
            << molecule.xmax*10.0 << " "
            << molecule.ymax*10.0 << " "
            << molecule.zmax*10.0 << std::endl;
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
        out << "#include \"toppar/" << molecule.name << ".itp\"" << std::endl;
    }
    out << std::endl;

    out << "[ system ]" << std::endl;
    out << "; Name" << std::endl;
    out << outputName << std::endl;
    out << std::endl;

    out << "[ molecules ]" << std::endl;
    out << "; Compound     #mols" << std::endl;
    for (auto& molecule : molecules)
    {
        out << molecule.name << "      " << molecule.count << std::endl;
    }
    out << std::endl;    

    out.close();
    
    return 0;
}