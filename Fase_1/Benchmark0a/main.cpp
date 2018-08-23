
#include "pzgmesh.h"
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>

#include "Monofasico.h"
#include "MonofasicoElastico.h"

#include "pzlog.h"
//------------------Benchmarks for Geomec------------------------


using namespace std;


int main(int argc, char *argv[])
{
    InitializePZLOG();
    MonofasicoElastico cenario0a;
    cenario0a.Run(2);

    return 0;
}



