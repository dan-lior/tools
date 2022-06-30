#include <iostream>
#include "affinity.h"

int main()
{
    AffinityIso<3> a;

    std::vector<double> xs({4.0, 6.0});

    std::cout <<  MiscellaneousMath::mean(xs) << std::endl;
}