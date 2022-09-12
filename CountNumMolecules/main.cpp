#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

std::map<char, float> masses;

void readMassesDB(string filename)
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
    return 0.0f;
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
            char element = line.c_str()[0];
            if (element != '#')
            {
            //std::string::size_type pos = line.find_first_of("1234567890");
            //std::string element = line.substr(0, pos);
                std::cout << "Element:" << element << std::endl;
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