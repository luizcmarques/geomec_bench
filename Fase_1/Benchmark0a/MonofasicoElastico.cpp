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
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "TPZInterfaceEl.h"
#include "TPZMatElastoPlasticDFN.h"
#include "TPZMatElastoPlasticDFN2D.h"
#include "TPZElasticCriterion.h"
#include "pzporoelastoplasticmem.h"
#include "pzcompelwithmem.h"
#include "TPZStiffFracture.h"
#include "pzelastoplasticanalysis.h"
#include "TPZDarcy2DMaterialMem.h"

#include "pzelasmat.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "pzbuildmultiphysicsmesh.h"
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
    fmatBCright=3;
    fmatBCtop=4;
    fmatBCleft=5;
    
    //Número de fraturas do problema:
    fnFrac = 1;
    
    fmatFrac.resize(fnFrac);
    fmatPointLeft.resize(fnFrac);
    fmatPointRight.resize(fnFrac);
    
    for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
        fmatFrac[i_frac] = 6+i_frac;
     //   fmatPointLeft[i_frac] = 7+3*i_frac;
     //   fmatPointRight[i_frac] = 8+3*i_frac;
    }
    
    //Material do elemento de interface
    fmatInterfaceLeft = 501;
    fmatInterfaceRight = 502;
    fmatFluxWrap= 503;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    
    //Materia de um ponto
    fmatPoint=-5;
    
    //Condições de contorno do problema
    fdirichlet =0;
    fneumann = 1;
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
    fEyoung = 1.;
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
    REAL Eyoung = 1.0;
    REAL poisson = 0.0;
    
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx = 0.;
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

    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(); //Função para criar a malha geometrica
   
    int n_div = 0;
    
    UniformRef(gmesh,n_div);
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    int p = 2;
    TPZCompEl::SetgOrder(p);
    
    // Malha Elásticidade
    TPZCompMesh *cmesh_E = CMesh_E(gmesh, pOrder); //Função para criar a malha computacional da velocidade
    
    // Malha Multifísica : Fluxo + Pressão
    
    TPZCompMesh *cmesh_q = CMesh_q(gmesh, pOrder); //Malha computacional de fluxo - Hdiv
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //Malha computacional de pressão - descontínua
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Malha multifísica (Fluxo+Pressão)
    
    
    
    std::ofstream filecEbf("MalhaC_E_before.txt"); //Impressão da malha computacional da velocidade (formato txt)
    cmesh_E->Print(filecEbf);
    
    std::vector<int> fracture_ids;
    for (int i_frac=0; i_frac< fnFrac; i_frac++) {
        fracture_ids.push_back(fmatFrac[i_frac]);
    }
            
    BreakH1Connectivity(*cmesh_E, fracture_ids); // Insert new connects to represent normal fluxes
    cmesh_E->ComputeNodElCon();
    {
        std::ofstream filecE("MalhaC_E.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_E->Print(filecE);
    }
    TPZManVector<TPZCompMesh *, 3> meshvector(1);
    meshvector[0] = cmesh_E;
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    std::ofstream filegvtk2("MalhaGeo2.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk2,true);

    //Resolução Analysis
    bool optimizeBandwidth = true;
    TPZAnalysis an(cmesh_E,optimizeBandwidth);
    TPZSymetricSpStructMatrix matskl(cmesh_E); //Pardiso
    //    TPZSkylineStructMatrix matskl(cmesh_E);
    int numthreads = 0;
    matskl.SetNumThreads(numthreads);
  //  TPZFMatrix<STATE> Initialsolution = an.Solution();
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
        std::ofstream filestiff("stiffness.txt");
//        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
//        an.Rhs().Print("Rhs = ",filerhs,EMathematicaInput);
#endif
    an.Solve();
    
//    {
//        //Pós-processamento standard (paraview):
//
//        std::string plotfile("Benchmark_0a_PoroElast_std.vtk");
//        TPZStack<std::string> scalnames, vecnames;
////        scalnames.Push("SigmaX");
////        scalnames.Push("SigmaY");
////        scalnames.Push("Pressure");
//        vecnames.Push("Displacement");
//
//        int postProcessResolution = 0; //  keep low as possible
//        int dim = gmesh->Dimension();
//        an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//        an.PostProcess(postProcessResolution,dim);
//    }
    
    
    // Acepting solution just because the problem is linear
    bool update_Q = true;
    { // for each material id
        int mat_id = 1;
        TPZMaterial *mat = an.Mesh()->FindMaterial(mat_id);
        if (!mat) {
            DebugStop();
        }
        
        // update mem
        {
            std::map<int, TPZMaterial *> & refMatVec = an.Mesh()->MaterialVec();
            std::map<int, TPZMaterial * >::iterator mit;
            
            TPZMatWithMem<TPZPoroElastoPlasticMem> * pMatWithMem2; // define in file pzporous.h
            
            for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
            {
                pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZPoroElastoPlasticMem> *>( mit->second);
                if(pMatWithMem2 != NULL)
                {
                    pMatWithMem2->SetUpdateMem(update_Q);
                }
            }
        }
        an.AssembleResidual();
    }
    
    TPZPostProcAnalysis post_an;
    post_an.SetCompMesh(an.Mesh());

    
    TPZManVector<int,10> post_mat_id(1);
    post_mat_id[0] = 1;
    TPZManVector<std::string,10> var_names(4), vecnames(0);
    var_names[0] = "XStress";
    var_names[1] = "YStress";
//    var_names[2] = "DisplacementMemX";
//    var_names[3] = "DisplacementMemY";
    var_names[2] = "DisplacementX";
    var_names[3] = "DisplacementY";
    post_an.SetPostProcessVariables(post_mat_id, var_names);
    
    {
        std::ofstream file("Malha_postproc.txt");
        post_an.Mesh()->Print(file);
    }
    
    TPZFStructMatrix structmatrix(post_an.Mesh());
    structmatrix.SetNumThreads(0);
    post_an.SetStructuralMatrix(structmatrix);
    post_an.TransferSolution();
    
    int dim = an.Mesh()->Dimension();
    int div = 0; //  keep low as possible
    std::string plotfile_post("Benchmark_0a_PoroElast.vtk");
    post_an.DefineGraphMesh(dim,var_names,vecnames,plotfile_post);
    post_an.PostProcess(div,dim);
    
    

//    int dim = gmesh->Dimension();
//    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//    an.PostProcess(postProcessResolution,dim);
//
//    std::cout << "FINISHED!" << std::endl;
//    {
//        std::ofstream filecE("MalhaC_E.txt"); //Impressão da malha computacional da velocidade (formato txt)
//        cmesh_E->Print(filecE);
//    }

}


void MonofasicoElastico::UniformRef(TPZGeoMesh * gmesh, int n_div){
    for ( int ref = 0; ref < n_div; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int n = gmesh->NElements();
        for ( int i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->Element(i);
            gel->Divide (filhos);
        }//for i
    }//ref
}

MonofasicoElastico::~MonofasicoElastico()
{
    
}

void MonofasicoElastico::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(1);
    const REAL Pi=M_PI;
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    STATE f_x = 0.;
    
    f[0] = f_x;
    
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
    
    grid = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/Benchmark0a/gmsh/GeometryBench.msh";

    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    Geometry.fPZMaterialId[0]["PointLeft"] = 0;
    Geometry.fPZMaterialId[0]["PointRight"] = 0;
    Geometry.fPZMaterialId[1]["bottom"] = fmatBCbott;
    Geometry.fPZMaterialId[1]["right"] = fmatBCright;
    Geometry.fPZMaterialId[1]["top"] = fmatBCtop;
    Geometry.fPZMaterialId[1]["left"] = fmatBCleft;
    Geometry.fPZMaterialId[1]["frac"] = fmatFrac[0];
    Geometry.fPZMaterialId[1]["frac2"] = fmatFrac[1];
    Geometry.fPZMaterialId[2]["Omega"] = fmatID;
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

    REAL y = ptx[1];
    sol.Resize(2, 0.);// ux, uy;
    
    sol[1] = 10.*y;
    deriv.Resize(2,2);//sigx, sigxy, sigyx, sigy
    deriv(0,0) = deriv(0,1) = deriv(1,0) = deriv(1,1) = 0.;
    
//    bool sol_dimensionless=true;
//
//    //REAL x = ptx[0];
//    REAL x = ptx[1];
//
//    REAL pini = 1000.;
//    REAL lamb = 8333.33;
//    REAL mi = 12500.0;
//    REAL H=1.;
//    REAL tp = ftimeatual;
//    int in;
//    REAL uD = 0.0, sigD=0.;
//    REAL sumuD = 0.0, sumsigD = 0.;
//
//    REAL M =0.;
//    REAL PI = atan(1.)*4.;
//
//    sol.Resize(2, 0.);// ux, uy;
//    deriv.Resize(2,2);//sigx, sigxy, sigyx, sigy
//    deriv(0,0) = deriv(0,1) = deriv(1,0) = deriv(1,1) = 0.;
//
//
//    REAL tD = tp;//(lamb+2.*mi)*perm*tp/(visc*H*H);
//    REAL xD = fabs(1.-x)/H;
//    for (in =999; in >= 0; in--) {
//
//        M = PI*(2.*in+1.)/2.;
//        sumuD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
//        sumsigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
//    }
//
//    uD = (H/H - xD) - sumuD;
//    sigD = -1. + sumsigD;
//
//    if(sol_dimensionless==true){
//        sol[1] = (-1.)*uD;
//        deriv(1,1) = sigD;
//    }else{
//        sol[1] = (-1.)*uD*(pini*H)/(lamb+2.*mi);
//        deriv(1,1) = (sigD)*pini;
//    }
    
}

TPZCompMesh *MonofasicoElastico::CMesh_E(TPZGeoMesh *gmesh, int pOrder)
{
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuousWithMem();
    
    //Material volumétrico
    int nstate = 2;
    TPZVec<REAL> force(fdim,0.);
    //REAL E = 0;
    //REAL poisson = 0;
    int planestress = -1;
    TPZMatElastoPlasticDFN2D<TPZElasticCriterion , TPZPoroElastoPlasticMem> *material;
    //material = new TPZElasticityMaterial(fmatID, fEyoung, fpoisson, ffx, ffy, planestress);
    material = new TPZMatElastoPlasticDFN2D<TPZElasticCriterion , TPZPoroElastoPlasticMem> (fmatID , planestress);
    TPZElasticCriterion obj ;
    
    TPZElasticResponse er;
    er.SetUp(fEyoung,fpoisson);
    obj.SetElasticResponse(er);
    material->SetPlasticity(obj);
    cmesh->InsertMaterialObject(material);
    
    //Material Fraturas
    TPZMat1dLin *materialFrac;
    for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
        materialFrac = new TPZMat1dLin(fmatFrac[i_frac]);
        TPZFMatrix<STATE> xkin3(nstate,nstate,0.), xcin3(nstate,nstate,0.), xbin3(nstate,nstate,0.), xfin3(nstate,nstate,0.);
        materialFrac->SetMaterial(xkin3, xcin3, xbin3, xfin3);
        cmesh->InsertMaterialObject(materialFrac);
    }
    
    //Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val1(1,1) = 1.e6;
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCbott, fmixed, val1, val2);
    val1.Zero();
    val2(1,0)= 1.;
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCtop, fneumann, val1, val2);
    val2.Zero();
    val2(1,0)= 0.;
    val1(0,0) = 1.e6;
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fmixed, val1, val2);
    val1.Zero();
    val2(0,0)=0.;
    TPZMaterial * BCond4 = material->CreateBC(material, fmatBCleft, fneumann, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
//    for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
//        //Mat Frac:
//        TPZMaterial * bc_frac_right = materialFrac->CreateBC(materialFrac, fmatPointRight[i_frac] , fdirichlet, val1, val2);
//        cmesh->InsertMaterialObject(bc_frac_right);
//        TPZMaterial * bc_frac_left = materialFrac->CreateBC(materialFrac, fmatPointLeft[i_frac] , fdirichlet, val1, val2);
//        cmesh->InsertMaterialObject(bc_frac_left);
//    }

    //Criando material para FluxWrap
//    TPZBndCond * bc_fracture_wrap;
//    bc_fracture_wrap = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
//    cmesh->InsertMaterialObject(bc_fracture_wrap);
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    // Inserir interface fratura - elem. volumétricos
    
  
    //Material Lagrange nas interfaces
    TPZStiffFracture *matInterLeft = new TPZStiffFracture(fmatInterfaceLeft, fdim, nstate);
    matInterLeft->SetMultiplier(-1);
    cmesh->InsertMaterialObject(matInterLeft);

    TPZStiffFracture *matInterRight = new TPZStiffFracture(fmatInterfaceRight, fdim, nstate);
    matInterRight->SetMultiplier(1);
    cmesh->InsertMaterialObject(matInterRight);
    
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
//    matids.insert(fmatFluxWrap);
    
//    cmesh->AutoBuild(matids);
//    gmesh->ResetReference();
//    matids.clear();
//    cmesh->LoadReferences();
//    cmesh->SetDefaultOrder(pOrder-1);
//    for (int i_frac = 0; i_frac < fnFrac; i_frac++) {
//        matids.insert(fmatFrac[i_frac]);
//        matids.insert(fmatPointLeft[i_frac]);
//        matids.insert(fmatPointRight[i_frac]);
//    }
//    matids.insert(fmatInterfaceRight);
//    matids.insert(fmatInterfaceLeft);
    
    cmesh->AutoBuild(matids);
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->InitializeBlock();
    
    
    return cmesh;
}

TPZCompMesh *MonofasicoElastico::CMesh_q(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
    
//    TPZDarcy2DMaterialMem<TPZPoroElastoPlasticMem> *material;
//    material = new TPZDarcy2DMaterialMem<TPZPoroElastoPlasticMem> (fmatID,fdim,1,1);
    
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
        
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
}


TPZCompMesh *MonofasicoElastico::CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    // @omar::
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Criando material:
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(fmatID);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    return cmesh;

}

TPZCompMesh *MonofasicoElastico::CMesh_m(TPZGeoMesh *gmesh, int pOrder){

    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    /// criar materiais
    
    //    TPZDarcy2DMaterialMem<TPZPoroElastoPlasticMem> *material;
    //    material = new TPZDarcy2DMaterialMem<TPZPoroElastoPlasticMem> (fmatID,fdim,1,1);
    
    
    TPZDarcy2DMaterialMem<TPZPoroElastoPlasticMem> *material = new TPZDarcy2DMaterialMem<TPZPoroElastoPlasticMem> (fmatID,fdim,1,1);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source);
//    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact);
    
    material->SetForcingFunction(fp);
   // material->SetForcingFunctionExact(solp);
    cmesh->InsertMaterialObject(material);
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fneumann, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    BCond0->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fneumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    BCond1->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fneumann, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    BCond2->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fneumann, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
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
        if (gel->MaterialId()!=fmatFrac[0]) {
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
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    
    for (unsigned int i_f = 0; i_f <  fracture_ids.size(); i_f++) {
        TPZFractureNeighborData fracture(cmesh.Reference(),fracture_ids[i_f],boundaries_ids);
        fracture.OpenFracture(&cmesh); // (ok)
        fracture.SetDiscontinuosFrac(&cmesh); // (ok)
        fracture.SetInterfaces(&cmesh, fmatInterfaceLeft, fmatInterfaceRight);
        
  //      std::ofstream filecE("CmeshWithFrac.txt"); //Impressão da malha computacional da velocidade (formato txt)
  //      cmesh.Print(filecE);
    }
    int aka  = 0;
}

void MonofasicoElastico::AddMultiphysicsInterfaces(TPZCompMesh &cmesh)
{

    DebugStop();
//    TPZGeoMesh *gmesh = cmesh.Reference();
//    std::set<int> velmatid;
//    velmatid.insert(fmatFrac[0]);
//
//
//    int64_t nel = gmesh->NElements();
//    for (int64_t el=0; el<nel; el++) {
//        TPZGeoEl *gel = gmesh->Element(el);
//        int matid = gel->MaterialId();
//        if(velmatid.find(matid) == velmatid.end())
//        {
//            continue;
//        }
//        int nsides = gel->NSides();
//        TPZGeoElSide gelside(gel,nsides-1);
//        TPZGeoElSide neighbour = gelside.Neighbour();
//        while (neighbour != gelside) {
//
//            TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
//            LeftElIndices[0]=0;
//            RightElIndices[0]=1;
//
//            if (neighbour.Element()->Dimension() == 1 && neighbour.Element()->MaterialId() == fmatFluxWrap) { //oioioi IDFlux -> ID
//                // create an interface element
//                TPZCompElSide celside = gelside.Reference();
//                TPZCompElSide Wrapneigh = neighbour.Reference();
//                if (!celside || !Wrapneigh) {
//                    DebugStop();
//                }
//                std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
//                " side " << neighbour.Side() <<
//                " and interface element " << gelside.Element()->Index() << std::endl;
//                TPZGeoElBC gelbc(gelside,fmatInterface);
//                int64_t index;
//                TPZMultiphysicsInterfaceElement *intf = new
//                TPZMultiphysicsInterfaceElement(cmesh,gelbc.CreatedElement(),index,Wrapneigh,celside);
//                intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
//
//            }
//            neighbour = neighbour.Neighbour();
//        }
//    }
    
}




