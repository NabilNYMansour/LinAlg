#include <iostream>
using namespace std;
#include "matrix.h"



int main(int argc, char const *argv[])
{
    matrix test(5);
    // test.makeIdentity();
    test.makeRandom(0,9);
    test.print();
    return 0;
}

