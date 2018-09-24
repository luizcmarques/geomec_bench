
#include "pzgmesh.h"
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>

#include "HidraulicoMonofasico2D.h"
#include "MonofasicoElastico.h"

#include "pzlog.h"
//------------------Benchmarks for Geomec------------------------

using namespace std;

/// Main dedicated to execute each scenario for benchmark_1
int main(int argc, char *argv[])
{
    InitializePZLOG();
    MonofasicoElastico scenario0a;
    scenario0a.Run(1);
    
    
//    MonofasicoElastico cenario0a;
//    cenario0a.Run(2);
    return 0;
}
