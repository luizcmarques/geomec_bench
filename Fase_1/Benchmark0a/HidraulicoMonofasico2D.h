/*
 *  HidraulicoMonofasico2D.h
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__HidraulicoMonofasico2D__
#define __PZ__HidraulicoMonofasico2D__

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "HidraulicoMonofasico2D.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "../TPZFractureInsertion.h"
#include "../TPZDarcyAnalysis.h"

using namespace std;
using namespace pzshape;

class HidraulicoMonofasico2D {
    
private:
    
    int fdim; //Dimensão do problema
    int fmatID; //Materia do elemento volumétrico
    int fdimFrac; //Dimensão da fratura
    
    //Materiais das condições de contorno
    int fmatBCbott;
    int fmatBCtop;
    int fmatBCleft;
    int fmatBCright;
    int fmatFrac;
    int fmatPointLeft;
    int fmatPointRight;
    
    //Material do elemento de interface
    int fmatInterface;
    int fmatFluxWrap;
    
    //Materiais das condições de contorno (elementos de interface)
    int fmatIntBCbott;
    int fmatIntBCtop;
    int fmatIntBCleft;
    int fmatIntBCright;
    
    //Materia de um ponto
    int fmatPoint;
    
    //Condições de contorno do problema
    int fdirichlet;
    int fneumann;
    int fpenetration;
    int fpointtype;
    int fdirichletPress;
    
    STATE fvisco;
    STATE fperm;
    int ftheta;
    
    bool finsert_fractures_Q = true;
    
public:

    HidraulicoMonofasico2D();

    void Run(int pOrder);
    
    ~HidraulicoMonofasico2D();
    
    /*  Malhas geometricas */
    TPZGeoMesh *CreateGMesh();

    /* Malhas computacionais */
    TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder);
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder);
    TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);
    
    //solucao exata
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //Fractures structure
    void Plot_over_fractures(TPZCompMesh *cmesh);
    void BreakConnectivity(TPZCompMesh &cmesh, int matId = 1);
    
    //Multiphysics Interfaces

    void AddMultiphysicsInterfaces(TPZCompMesh &cmesh);
    
};


#endif 
