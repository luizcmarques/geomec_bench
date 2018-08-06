/*
 *  Monofasico.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "MonofasicoElastico.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzelasmat.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"

#include "pzporoelasticmf2d.h"

#include "TPZDarcy1DMaterial.h"
#include "TPZDarcy2DMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"


#define TRIANGLEMESH

using namespace std;

const REAL Pi=M_PI;
REAL ftimeatual = 0.;


MonofasicoElastico::MonofasicoElastico()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    fdimFrac = 1; //Dimensão da fratura
    
    //Materiais das condições de contorno
    fmatBCbott=2;
    fmatBCtop=3;
    fmatBCright=4;
    fmatBCleft=5;
    fmatFrac=6;
    fmatPointLeft = 7;
    fmatPointRight = 8;
    
    //Material do elemento de interface
    fmatInterface=17;
    fmatFluxWrap=21;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    
    //Materia de um ponto
    fmatPoint=-5;
    
    //Condições de contorno do problema
    fdirichlet =0;
    fneumann = 11;
    fpenetration=2;
    fpointtype=5;
    fdirichletPress=6;
    fneumdir=10;
    fdirfreey_neum=300;
    fdirneum = 1;
    fmixedneum = 21;
    fmixeddirich = 20;
    fmixedFreeYXdirich = 400;
    fmixedFreeXYdirich = 500;
    
    ftheta=-1;
    fEyoung = 0.;
    fpoisson = 0.;
    falpha = 0.;
    fSe = 0.;
    fperm = 0.;
    fvisc = 0.;
    ffx = 0.;
    ffy = 0.;
    fsign = 0.;
    
    fpref = 0.;
    fLref = 0.;
    fkovervisc = 0.;
    fvalsourceterm = 0.;
    
}


void MonofasicoElastico::Run(int pOrder)
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    HDivPiola = 1;
    
    TPZMaterial::gBigNumber = 1.e16;
    REAL Eyoung = 3.e6;
    REAL poisson = 0.2;
    
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx = 0.0;
    REAL fy = gravity*rockrho;
    
    REAL alpha = 0.;
    REAL Se = 0.0;
    REAL perm = 1.;
    REAL visc = 1.;
    REAL sig0 = 0.;
    REAL pini = 0.;
    REAL La = 1.;
    REAL Ly = 5.*La;
    REAL Lx = 8.*La;
    REAL timeT = 0.;
    this->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
    
    ofstream saidaerro("ErroLoula.txt");
    
    for(int p = 2; p < 3; p++)
    {
        int pu = p;
        int pq = pu;
        int pp;
//        if(triang==true){
//            pp = pq-1;
//        }else{
            pq=pu-1;
            pp = pq;
//        }
        
        int h;
        saidaerro<<"\n CALCULO DO ERRO, ELEM. TRIANG., COM ORDEM POLINOMIAL pu = "<< pu << ", pq = "<< pq << " e pp = "<< pp<<endl;
        for (h = 3; h< 4; h++)
        {
            
            saidaerro<<"\n========= PARA h = "<< h<<"  ============= "<<endl;
    

    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:

    TPZCompMesh *cmesh_E = CMesh_E(gmesh, pOrder); //Função para criar a malha computacional da velocidade
  //  TPZCompMesh *cmesh_q = CMesh_q(gmesh, pOrder); //Função para criar a malha computacional da velocidade
  //  TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //Função para criar a malha computacional da pressão
    
    {
        std::ofstream filecE("MalhaC_E.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecq("MalhaC_q.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_E->Print(filecE);
 //       cmesh_q->Print(filecq);
 //       cmesh_p->Print(filecp);
    }
    
    std::vector<int> fracture_ids;
    fracture_ids.push_back(fmatFrac);
    BreakH1Connectivity(*cmesh_E, fracture_ids); // Insert new connects to represent normal fluxes
    TPZPoroElasticMF2d *mymaterial;
    
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder, mymaterial); //Função para criar a malha computacional multifísica
    
    TPZManVector<TPZCompMesh *, 3> meshvector(1);
    meshvector[0] = cmesh_E;
  //  meshvector[1] = cmesh_q;
  //  meshvector[2] = cmesh_p;
    
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    {
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
    
//    AddMultiphysicsInterfaces(*cmesh_m);
    
    {
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
    
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    //Resolvendo o Sistema:
  //  int NDeltaT = 1000000;
  //  int intervsaidas = NDeltaT/20;
 //   REAL deltaT=timeT/NDeltaT; //second
 //   mymaterial->SetTimeStep(deltaT);
 //   REAL maxTime = timeT;
    
    TPZAnalysis an(cmesh_m);
    TPZFMatrix<STATE> Initialsolution = an.Solution();
    //Initialsolution.Print("solini");
    
    std::string outputfile;
    outputfile = "TransientSolution";
    int nthreads = 0;
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = this->MassMatrix(mymaterial, cmesh_m, nthreads);
    
    //Criando matriz de rigidez (matK) e vetor de carga
    TPZFMatrix<STATE> matK;
    TPZFMatrix<STATE> fvec;
    this->StiffMatrixLoadVec(mymaterial, cmesh_m, an, matK, fvec, nthreads);
    
    int nrows;
    nrows = matM->Rows();
    TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
    TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
    TPZFMatrix<STATE> Lastsolution = Initialsolution;
    

        an.Solve();

            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            this->PosProcessMultphysics(meshvector,cmesh_m,an,plotfile);

            TPZVec<REAL> erros;
            
            saidaerro<<" Erro da simulacao multifisica do deslocamento (u)" <<endl;
            TPZAnalysis an12(cmesh_E);
            an12.SetExact(*Sol_exact);
       //     an12.PostProcessError(erros, false, saidaerro);
            
            saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
       //     TPZAnalysis an22(cmesh_q);
       //     an22.SetExact(*Sol_exact);
       //     an22.PostProcessError(erros, false, saidaerro);
            
            saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
        //    TPZAnalysis an32(cmesh_p,false);
        //    an32.SetExact(*Sol_exact);
        //    an32.PostProcessError(erros, false, saidaerro);


    
    cmesh_E->CleanUp();
 //   cmesh_q->CleanUp();
 //   cmesh_p->CleanUp();
    //mphysics->CleanUp();
    delete cmesh_E;
 //   delete cmesh_q;
 //   delete cmesh_p;
    //delete mphysics;
    delete gmesh;
    }
}

//    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
//    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
//    TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso nao simetrico ***
//    matskl.SetNumThreads(numthreads);
//
//    an.SetStructuralMatrix(matskl);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELU);
//    an.SetSolver(step);
//    an.Assemble(); //Assembla a matriz de rigidez (e o vetor de carga) global
//    an.Solve();
//
//    //Pós-processamento (paraview):
//
//    std::string plotfile("Benchmark_0a_DarcyTest.vtk");
//    TPZStack<std::string> scalnames, vecnames;
//    scalnames.Push("P");
//    vecnames.Push("V");
//
//    int postProcessResolution = 0; //  keep low as possible
//    int dim = gmesh->Dimension();
//    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//    an.PostProcess(postProcessResolution,dim);
//
//    if(insert_fractures_Q) {
//        std::cout << "Postprocessing fracture" << std::endl;
//        Plot_over_fractures(cmesh_m);
//    }
//
    std::cout << "FINISHED!" << std::endl;
    
}


MonofasicoElastico::~MonofasicoElastico()
{
    
}

TPZAutoPointer <TPZMatrix<STATE> > MonofasicoElastico::MassMatrix(TPZPoroElasticMF2d * mymaterial, TPZCompMesh* mphysics, int nthreads){
    
    mymaterial->SetLastState();
    //TPZSkylineStructMatrix matsp(mphysics);
    //TPZSkylineNSymStructMatrix matsp(mphysics);
    TPZSpStructMatrix matsp(mphysics);
    matsp.SetNumThreads(nthreads);
    
    std::set< int > materialid;
    int matid = mymaterial->MatId();
    materialid.insert(matid);
    matsp.SetMaterialIds (materialid);
    TPZAutoPointer<TPZGuiInterface> guiInterface;
    TPZFMatrix<STATE> Un;
    
    
    //TPZMatrix<REAL> *matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

void MonofasicoElastico::StiffMatrixLoadVec(TPZPoroElasticMF2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<STATE> &matK1, TPZFMatrix<STATE> &fvec, int nthreads) {
    
    mymaterial->SetCurrentState();
    //TPZFStructMatrix matsk(mphysics);
    TPZSkylineStructMatrix matsk(mphysics);
    matsk.SetNumThreads(nthreads);
    an.SetStructuralMatrix(matsk);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
    an.SetSolver(step);
    an.Run();
    
    matK1 = an.StructMatrix();
    fvec = an.Rhs();
    
}

void MonofasicoElastico::SetParameters(REAL mod_young, REAL mod_poisson, REAL coef_alpha, REAL coef_Se, REAL permeabil_fluido, REAL visc_fluido, REAL fx, REAL fy,REAL sign){
    
    fEyoung = mod_young;
    fpoisson= mod_poisson;
    falpha = coef_alpha;
    fSe = coef_Se;
    fperm = permeabil_fluido;
    fvisc = visc_fluido;
    ffx = fx;
    ffy = fy;
    fsign = sign;
}

void MonofasicoElastico::PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    
    TPZManVector<std::string,10> scalnames(3), vecnames(2);
    //scalnames[0] = "DisplacementX";
    //scalnames[1] = "DisplacementY";
    vecnames[0] = "Displacement";
    scalnames[0] = "SigmaX";
    scalnames[1] = "SigmaY";
    //scalnames[0] = "PorePressure";
    //scalnames[2] = "FluxoY";
    //vecnames[0] = "Fluxo";
    //vecnames[1] = "MinusKMuGradP";
    
    scalnames[2] = "ExactPressure";
    //scalnames[6] = "FluxoX";
    
    //scalnames[4] = "ExactDisplacementX";
    //scalnames[5] = "ExactDisplacementY";
    // scalnames[8] = "ExactSigmaX";
    //scalnames[6] = "ExactSigmaY";
    vecnames[1]  = "ExactFluxo";
    //vecnames[1]  = "ExactDisplacement";
    //vecnames[4] = "MinusKMuGradP";
    
    const int dim = 2;
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    //    std::ofstream out("malha.txt");
    //    an.Print("nothing",out);
}

TPZGeoMesh *MonofasicoElastico::CreateGMesh()
{
    int64_t id, index;
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Aqui é implementado um método para malhas criadas no GMSH
    
    //std::string dirname = PZSOURCEDIR;
    std::string grid;
    
    grid = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/Benchmark0a/gmsh/GeometryBench0b.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    gmesh = Geometry.GeometricGmshMesh(grid);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //    int n_div = 0;
    //    UniformRefine(gmesh,n_div);
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

}

void MonofasicoElastico::Sol_exact(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv){
    
    bool sol_dimensionless=true;
    
    //REAL x = ptx[0];
    REAL x = ptx[1];
    
    REAL pini = 1000.;
    REAL lamb = 8333.33;
    REAL mi = 12500.0;
    REAL H=1.;
    REAL tp = ftimeatual;
    int in;
    REAL uD = 0.0, sigD=0.;
    REAL sumuD = 0.0, sumsigD = 0.;
    
    REAL M =0.;
    REAL PI = atan(1.)*4.;
    
    sol.Resize(2, 0.);// ux, uy;
    deriv.Resize(2,2);//sigx, sigxy, sigyx, sigy
    deriv(0,0) = deriv(0,1) = deriv(1,0) = deriv(1,1) = 0.;
    
    
    REAL tD = tp;//(lamb+2.*mi)*perm*tp/(visc*H*H);
    REAL xD = fabs(1.-x)/H;
    for (in =999; in >= 0; in--) {
        
        M = PI*(2.*in+1.)/2.;
        sumuD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
        sumsigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
    }
    
    uD = (H/H - xD) - sumuD;
    sigD = -1. + sumsigD;
    
    if(sol_dimensionless==true){
        sol[1] = (-1.)*uD;
        deriv(1,1) = sigD;
    }else{
        sol[1] = (-1.)*uD*(pini*H)/(lamb+2.*mi);
        deriv(1,1) = (sigD)*pini;
    }
    
}

TPZCompMesh *MonofasicoElastico::CMesh_E(TPZGeoMesh *gmesh, int pOrder)
{
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsContinuous();
    
    /// criar material
    TPZVec<REAL> force(fdim,0.);
    //REAL E = 0;
    //REAL poisson = 0;
    int planestress = -1;
    TPZElasticityMaterial *material;
    
    material = new TPZElasticityMaterial(fmatID, fEyoung, fpoisson, ffx, ffy, planestress);
    cmesh->InsertMaterialObject(material);
    
    TPZMat1dLin *materialFrac = new TPZMat1dLin(fmatFrac);
    cmesh->InsertMaterialObject(materialFrac);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    //Mat Frac:
    
    TPZMaterial * bc_frac_right = materialFrac->CreateBC(materialFrac, fmatPointRight , fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(bc_frac_right); //Insere material na malha
    TPZMaterial * bc_frac_left = materialFrac->CreateBC(materialFrac, fmatPointLeft , fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(bc_frac_left); //Insere material na malha
    //Criando material para FluxWrap
    
    TPZBndCond * bc_fracture_wrap;
    bc_fracture_wrap = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
    cmesh->InsertMaterialObject(bc_fracture_wrap);
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    std::set<int> matids;
    matids.insert(fmatID);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCright);
    matids.insert(fmatBCtop);
    matids.insert(fmatBCleft);
    matids.insert(fmatFluxWrap);
    
    cmesh->AutoBuild(matids);
    gmesh->ResetReference();
    matids.clear();
    matids.insert(fmatFrac);
    matids.insert(fmatPointLeft);
    matids.insert(fmatPointRight);
    
    if (insert_fractures_Q) {
        cmesh->SetDimModel(fdimFrac);
        cmesh->SetAllCreateFunctionsHDiv();
        cmesh->AutoBuild(matids);
        cmesh->AdjustBoundaryElements();
        cmesh->ExpandSolution();
    }
    return cmesh;
}

TPZCompMesh *MonofasicoElastico::CMesh_q(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(fmatID,fdim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2);
    
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}


TPZCompMesh *MonofasicoElastico::CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    
    /// criar materiais
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(fmatID,fdim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fdim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();

    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    ///inserir connect da pressao
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
//        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//        {
//            if(triang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
//            else celdisc->SetTensorialShape();
//        }
    }
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;

}

TPZCompMesh *MonofasicoElastico::CMesh_m(TPZGeoMesh *gmesh, int pOrder, TPZPoroElasticMF2d * &mymaterial){

    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=fmatID;
    
    //    criar material
    int planestress = 1; // This is a Plain strain problem
    
    mymaterial = new TPZPoroElasticMF2d(MatId,fdim);
    mymaterial->SetfPlaneProblem(planestress);
    
    mymaterial->SetParameters(fperm, fvisc);
    mymaterial->SetElasticityParameters(fEyoung, fpoisson, falpha, fSe, ffx, ffy);
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(Sol_exact);
    mymaterial->SetForcingFunctionExact(solExata);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
//    REAL sig0 = fsign;
    REAL sig0 = 0.;
    REAL ptop = 0.;
    val2(0,0)= 0.;
    val2(1,0)= 1.;
    val2(2,0)= 0.;
    TPZMaterial * BCond1 = mymaterial->CreateBC(mat, fmatBCbott,fneumann, val1, val2);
//
    val2(0,0)= 0.;
    val2(1,0)= 0.;
    val2(2,0)= 0.;
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat,fmatBCtop,fdirichlet, val1, val2);

 
    val2(0,0)= 1.;
    val2(1,0)= 0.;
    val2(2,0)= 0.;
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat,fmatBCright, fneumann, val1, val2);
    
    val2.Redim(3,1);
    val2(0,0)= 0.;
    val2(1,0)=0.;
    REAL big = mymaterial->gBigNumber;
    val1(0,0) = 0.;
 
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat,fmatBCleft, fdirichlet, val1, val2);
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    
    //Para ser criado depois (main)
    // Creating multiphysic elements into mphysics computational mesh
//    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
//    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);

    return mphysics;
    
}


void MonofasicoElastico::Plot_over_fractures(TPZCompMesh *cmesh){
    // Counting for n_fractures
    int ncel = cmesh->NElements();
    TPZStack<int> frac_indexes;
    for (int icel = 0; icel<ncel; icel++) {
        TPZCompEl *cel=cmesh->Element(icel);
        if (!cel) {
            DebugStop();
        }
        TPZGeoEl *gel=cel->Reference();
        if (!gel) {
            DebugStop();
        }
        if (gel->MaterialId()!=fmatFrac) {
            continue;
        }
        frac_indexes.Push(cel->Index());
    }
    
    int n_frac_cels = frac_indexes.size();
    int var_p = 0;
    TPZVec<REAL> par_xi(1,0.0);
    TPZManVector<STATE,3> sol,x(3,0.0);
    TPZFMatrix<REAL> pressure(n_frac_cels*3,3,0.0);
    
    TPZManVector<REAL,3> par_vals(3,0.0);
    par_vals[0] = -1.0;
    par_vals[1] =  0.0;
    par_vals[2] = +1.0;
    
    
    for (int ifrac = 0; ifrac < n_frac_cels; ifrac++) {
        
        TPZCompEl *cel= cmesh->Element(frac_indexes[ifrac]);
        TPZGeoEl *gel=cel->Reference();
        for (int ip = 0; ip < par_vals.size(); ip++) {
            par_xi[0] = par_vals[ip];
            gel->X(par_xi, x);
            cel->Solution(par_xi, var_p, sol);
            pressure(ifrac*3+ip,0) = x[0];
            pressure(ifrac*3+ip,1) = x[1];
            pressure(ifrac*3+ip,2) = sol[0];
        }
        
    }
    
    pressure.Print("pf = ",std::cout,EMathematicaInput);
}


void MonofasicoElastico::BreakConnectivity(TPZCompMesh &cmesh, int matId)
{
    
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    cmesh.LoadReferences();
    cmesh.SetDimModel(fdim);
    cmesh.SetAllCreateFunctionsHDiv();
    int64_t ncel = cmesh.NElements();
    
    
    for (int64_t el=0; el<ncel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId()!=matId){
            continue;
        }
        if(gel->Dimension()!=fdim-1){
            DebugStop();
        }
        
        TPZStack<TPZCompElSide> neigh;
        int nsides = gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        
        gelside.EqualLevelCompElementList(neigh, 0, 0);
        
        if(neigh.size()!=2){
            DebugStop();
        }
        gel->ResetReference();
        neigh[0].Element()->Reference()->ResetReference();
        neigh[1].Element()->Reference()->ResetReference();
        
        //working on element 0
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neigh[0].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            
            int locindex = intel->MidSideConnectLocId(neigh[0].Side());
            TPZConnect &midsideconnect = intel->MidSideConnect(neigh[0].Side());
            if(midsideconnect.NElConnected() != 2)
            {
                DebugStop();
            }
            
            //Duplica um connect
            int64_t index = cmesh.AllocateNewConnect(midsideconnect.NShape(), midsideconnect.NState(), midsideconnect.Order());
            
            intel->SetConnectIndex(locindex, index);
            midsideconnect.DecrementElConnected();
            cmesh.ConnectVec()[index].IncrementElConnected();
            intel->SetSideOrient(neigh[0].Side(), 1);
            
            
            
            TPZGeoElBC bc(intel->Reference(),neigh[0].Side(),fmatFluxWrap);
            cmesh.CreateCompEl(bc.CreatedElement(), index);
            
            TPZCompEl *var = cmesh.Element(index);
            var->Reference()->ResetReference();
            intel->Reference()->ResetReference();
            
            
        }
        
        // working on element 1
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neigh[1].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            
            intel->SetSideOrient(neigh[1].Side(), 1);
            
            int64_t index;
            
            TPZGeoElBC bc(intel->Reference(),neigh[1].Side(),fmatFluxWrap);
            cmesh.CreateCompEl(bc.CreatedElement(), index);
            TPZCompEl *var = cmesh.Element(index);
            var->Reference()->ResetReference();
            
            intel->Reference()->ResetReference();
            
        }
        
    }
    
    cmesh.ExpandSolution();
}

void MonofasicoElastico::BreakH1Connectivity(TPZCompMesh &cmesh, std::vector<int> fracture_ids)
{
    for (unsigned int i_f = 0; i_f <  fracture_ids.size(); i_f++) {
        TPZFractureNeighborData fracture(cmesh.Reference(),fracture_ids[i_f]);
        int aka = 0;
    }
}


void MonofasicoElastico::AddMultiphysicsInterfaces(TPZCompMesh &cmesh)
{
    
    TPZGeoMesh *gmesh = cmesh.Reference();
    std::set<int> velmatid;
    velmatid.insert(fmatFrac);
    
    
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if(velmatid.find(matid) == velmatid.end())
        {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            
            TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
            LeftElIndices[0]=0;
            RightElIndices[0]=1;
            
            if (neighbour.Element()->Dimension() == 1 && neighbour.Element()->MaterialId() == fmatFluxWrap) { //oioioi IDFlux -> ID
                // create an interface element
                TPZCompElSide celside = gelside.Reference();
                TPZCompElSide Wrapneigh = neighbour.Reference();
                if (!celside || !Wrapneigh) {
                    DebugStop();
                }
                std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                " side " << neighbour.Side() <<
                " and interface element " << gelside.Element()->Index() << std::endl;
                TPZGeoElBC gelbc(gelside,fmatInterface);
                int64_t index;
                TPZMultiphysicsInterfaceElement *intf = new
                TPZMultiphysicsInterfaceElement(cmesh,gelbc.CreatedElement(),index,Wrapneigh,celside);
                intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                
            }
            neighbour = neighbour.Neighbour();
        }
    }
    
}




