/*
 *  Monofasico.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "Monofasico.h"
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

#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"

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


Monofasico::Monofasico()
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
    fdirichlet=0;
    fneumann=1;
    fpenetration=2;
    fpointtype=5;
    fdirichletPress=6;
    
    fvisco=1;
    fperm=1;
    ftheta=-1;
    
}


void Monofasico::Run(int pOrder)
{

    TPZMaterial::gBigNumber = 1.e16;
    
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    HDivPiola = 1;
    
    pOrder = 2; //Ordem polinomial de aproximação
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(); //Função para criar a malha geometrica
    
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_v = CMesh_v(gmesh, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //Função para criar a malha computacional da pressão
    
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
    }
    
    BreakConnectivity(*cmesh_v, fmatFrac); // Insert new connects to represent normal fluxes
    
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_v->Print(filecv);
    }
    
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Função para criar a malha computacional multifísica
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    {
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
    
    AddMultiphysicsInterfaces(*cmesh_m);
    
    {
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
    
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    matskl.SetNumThreads(numthreads);
    
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.Assemble(); //Assembla a matriz de rigidez (e o vetor de carga) global
    an.Solve();
    
    //Pós-processamento (paraview):
    
    std::string plotfile("Benchmark_0a_DarcyTest.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    
    int postProcessResolution = 0; //  keep low as possible
    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    if(insert_fractures_Q) {
        std::cout << "Postprocessing fracture" << std::endl;
        Plot_over_fractures(cmesh_m);
    }
    
    std::cout << "FINISHED!" << std::endl;
    
}


Monofasico::~Monofasico()
{
    
}

TPZGeoMesh *Monofasico::CreateGMesh()
{
    int64_t id, index;
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Aqui é implementado um método para malhas criadas no GMSH
    
    //std::string dirname = PZSOURCEDIR;
    std::string grid;
    
    grid = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/Benchmark0a/gmsh/msh/GeometryBench.msh";
    
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

void Monofasico::Sol_exact(const TPZVec<REAL> & x, TPZVec<STATE>& sol, TPZFMatrix<STATE> &dsol){
    
    dsol.Resize(3,2);
    sol.Resize(3);
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    sol.resize(3);
    
    STATE theta = atan(yv/xv);
    STATE r=sqrt(xv*xv+yv*yv);
    
    STATE v_x =  -r*sin(theta);
    STATE v_y =  r*cos(theta);
    STATE p =   0.;
    
    sol[0] = v_x; // x direction
    sol[1] = v_y; // y direction
    sol[2] = p; //
    
}

TPZCompMesh *Monofasico::CMesh_v(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID); //Criando material que implementa a formulação fraca do problema modelo
    TPZMat1dLin *materialFrac = new TPZMat1dLin(fmatFrac);
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->InsertMaterialObject(materialFrac); //Insere material na malha
    
    //cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
    //cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xbin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    materialFrac->SetMaterial(xkin, xcin, xbin, xfin);
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2);//Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //Mat Frac:
    
    TPZMaterial * BCond4 = materialFrac->CreateBC(materialFrac, fmatPointRight , fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    
    TPZMaterial * BCond5 = materialFrac->CreateBC(materialFrac, fmatPointLeft , fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha
    
    //Criando material para FluxWrap
    
    TPZBndCond *FluxWrapBC;
    FluxWrapBC = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
    cmesh->InsertMaterialObject(FluxWrapBC);
    
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
    
    cmesh->SetDimModel(fdim);
    return cmesh;
    
}


TPZCompMesh *Monofasico::CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //pOrder--; // Space restriction apapapa
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    
    //Definição do espaço de aprximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID);//criando material que implementa a formulacao fraca do problema modelo
    TPZMat2dLin *materialFrac = new TPZMat2dLin(fmatFrac);//criando material que implementa a formulacao fraca do problema modelo
    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->InsertMaterialObject(materialFrac); //Insere material na malha
    cmesh->SetAllCreateFunctionsDiscontinuous(); //Criando funções
    
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    materialFrac->SetMaterial(xkin, xcin, xfin);
    
    
    std::set<int> matids;
    matids.insert(fmatID);
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    cmesh->AutoBuild(matids);
    
    if (insert_fractures_Q) {
        gmesh->ResetReference();
        matids.clear();
        matids.insert(fmatFrac);
        cmesh->SetDimModel(fdimFrac); //Insere dimensão do modelo
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->AutoBuild(matids);
    }
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    cmesh->SetDimModel(fdim);
    return cmesh;
    
}

TPZCompMesh *Monofasico::CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    //TPZPoroElasticMF2d *material = new TPZPoroElasticMF2d(matID,dim);
    TPZDarcy2DMaterial *material = new TPZDarcy2DMaterial(fmatID,fdim,1,ftheta);//criando material que implementa a formulacao fraca do problema modelo
    TPZFMatrix<REAL> K(fdim,fdim),invK(fdim,fdim);
    K.Zero();
    invK.Zero();
    
    K(0,0)=3.38801e-7;
    K(1,1)=2.566e-11;
    invK(0,0)=1./K(0,0);
    invK(1,1)=1./K(1,1);
    
    // material->SetPermeabilityTensor(K, invK);
    
    TPZDarcy2DMaterial *materialFrac = new TPZDarcy2DMaterial(fmatFrac,fdimFrac,1,ftheta);//criando material que implementa a formulacao fraca do problema modelo
    REAL kf = 4.68789e-4;
    REAL Dyf = 6.5e-5;
    materialFrac->SetPermeability(kf*Dyf);
    
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact);
    
    material->SetForcingFunctionExact(solp);
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(materialFrac);
    
    
    //Condições de contorno:
    
    STATE Pjusante = 54.9;
    STATE Pmontante = 55.0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
    val1(0,0) = 0.; //botton
    
    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fneumann , val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(sol_exact1,bc_inte_order);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    val1(0,0) = 0.; //top
    
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fneumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact2,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    val1(0,0) = Pjusante; // right
    
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact1,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    val1(0,0) = Pmontante; // left
    
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond2->SetForcingFunction(solucao_exact2,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    if (insert_fractures_Q) {
        val1(0,0) =  -Pjusante; // right
        
        TPZMaterial * BCond4 = materialFrac->CreateBC(materialFrac, fmatPointRight , fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        cmesh->InsertMaterialObject(BCond4); //Insere material na malha
        
        val1(0,0) = -Pmontante; // left
        
        TPZMaterial * BCond5 = materialFrac->CreateBC(materialFrac, fmatPointLeft , fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        cmesh->InsertMaterialObject(BCond5); //Insere material na malha
    }
    
    val1(0,0) = 0.0;
    
    if (insert_fractures_Q) {
        TPZMaterial *MatLagrange = new TPZLagrangeMultiplier(fmatInterface,fdimFrac,1); //
        cmesh->InsertMaterialObject(MatLagrange);
        
        TPZBndCond *FluxWrapBC = material->CreateBC(material,fmatFluxWrap,fdirichlet,val1,val2);
        cmesh->InsertMaterialObject(FluxWrapBC);
    }
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

void Monofasico::Plot_over_fractures(TPZCompMesh *cmesh){
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


void Monofasico::BreakConnectivity(TPZCompMesh &cmesh, int matId)
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

void Monofasico::AddMultiphysicsInterfaces(TPZCompMesh &cmesh)
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




