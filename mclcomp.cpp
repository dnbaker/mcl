#include "mcl.h"

int main() {
    blaze::DynamicMatrix<double> mat = blaze::generate(10, 10, [](auto,auto) {return std::rand() * 0.002344;});
    bmcl::mcl(mat, bmcl::MCLSettings());
}
