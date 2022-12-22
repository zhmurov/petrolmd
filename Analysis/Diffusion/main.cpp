#include <iostream>
#include <chemfiles.hpp>

int main() {
    chemfiles::Trajectory file("filename.xyz");

    auto frame = file.read();
    std::cout << "There are " << frame.size() << " atoms in the frame" << std::endl;
    auto positions = frame.positions();

    // Do awesome science here with the positions

    if (frame.velocities()) {
        auto velocities = *frame.velocities();

        // If the file contains information about the
        // velocities, you will find them here.
    }

    return 0;
}