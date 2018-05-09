//
//  TPZFracSimulation.cpp
//  PZ
//
//  Created by Philippe Devloo on 26/09/17.
//

#include "TPZFracSimulation.h"

#include "mixedpoisson.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZMixedPoissonParabolic.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"

#include "TPZRefPattern.h"

#include <algorithm>
using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.fracsimulation"));
#endif


/// Build an MHM object with information given by the root name
void TPZFracSimulation::ReadDataFile(const std::string &rootname)
{
    std::string datafilename = rootname + ".data";
    std::ifstream datafile(datafilename);
    int dimension = fMHM->GMesh()->Dimension();
    ReadDataFile(datafile);
    std::string meshfilename = rootname + ".msh";
    fGmsh.GeometricGmshMesh(meshfilename,fMHM->GMesh().operator->());
#ifdef PZDEBUG
    if(fMHM->GMesh()->NElements() == 0)
    {
        std::cout << "Something went wrong with reading " << meshfilename << std::endl;
        DebugStop();
    }
#endif
    SetStandardBoundaryElements();

    AdjustGeometricMesh(rootname);
    // in order to build the wrap mesh
    fMHM->DivideSkeletonElements(0);
#ifdef PZDEBUG
    {
        std::string filename = rootname + "_gmesh_wrap.vtk";
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh(), out);
    }
#endif
}

/// read the data in the data file and put the information in right places
// creates material objects for the computational meshes
void TPZFracSimulation::ReadDataFile(std::ifstream &input)
{
    ReadPreamble(input);
    ReadFractures(input);
}

static void ReadNextLine(std::ifstream &input, std::string &line)
{
    std::getline(input,line);
    while (input && line[0] == '#') {
        std::getline(input,line);
    }
}

static std::string stripquotes(const std::string &input)
{
    char first = input.front();
    char last = input.back();
    int firstindex = 0;
    std::size_t length = input.size();
    if(first == '\"')
    {
        firstindex++;
        length--;
    }
    if(last == '\"') length--;
    string result = input.substr(firstindex,length);
    return result;
}

/// reads the preamble of the data file
void TPZFracSimulation::ReadPreamble(std::ifstream &input)
{

    int dimension = fMHM->GMesh()->Dimension();
    std::string line;
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        if (dimension == 2) {
            stin >> fFracSet.fLowLeft[0] >> fFracSet.fLowLeft[1] >> fFracSet.fTopRight[0] >> fFracSet.fTopRight[1];
        }
        else if(dimension == 3)
        {
            fFracSet.fLowLeft.Resize(3, 0.);
            fFracSet.fTopRight.Resize(3, 0.);
            stin >> fFracSet.fLowLeft[0] >> fFracSet.fLowLeft[1] >> fFracSet.fLowLeft[2];
            REAL tmp;
            for (int i=0; i<15; i++) {
                stin >> tmp;
            }
            stin >> fFracSet.fTopRight[0] >> fFracSet.fTopRight[1] >> fFracSet.fTopRight[2];
        }
    }
    REAL delxmin = min(fFracSet.fTopRight[0]-fFracSet.fLowLeft[0],fFracSet.fTopRight[1]-fFracSet.fLowLeft[1]);
    REAL tol = delxmin/200.; // will generate a 200x200 raster for points
    fFracSet.SetTol(tol);
    ReadNextLine(input, line);
    int numMHM;
    {    std::istringstream stin(line);
        stin >> numMHM;
    }
    fFracSet.SetMHMSpacing(numMHM);
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fFracSet.fElementSize >> fFracSet.fMinElementSize;
    }
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fSimulationType >> fInitialPressure;
        if(fSimulationType != 0 && fSimulationType != 1) DebugStop();
    }
    ReadNextLine(input, line);
    int nummat;
    int matidcounter = 20;
    {
        std::istringstream stin(line);
        stin >> nummat;
    }
    fGmsh.fPZMaterialId[1]["MHMLine"] = fMHM->fSkeletonMatId;
    fGmsh.fPZMaterialId[2]["MHMSurface"] = fMHM->fSkeletonMatId;
    std::string planemat;
    for (int i=0; i<nummat; i++) {
        ReadNextLine(input, line);
        std::string matname , mattype;
        int matid, matdimension;
        REAL density, perm;
        {
            std::istringstream stin(line);
            stin >> matname >> matdimension >> matid >> mattype >> density >> perm;
        }
        matname.erase(0,1);
        matname.erase(matname.end()-1,matname.end());
        if (mattype == "flow" && matdimension == dimension) {
            fFracSet.fPhysicalname = matname;
            fGmsh.fPZMaterialId[dimension][matname] = matidcounter;
            InsertDarcyMaterial(matidcounter, perm, density);
            fMaterialIds[matname] = matidcounter;
            matidcounter++;
        }
        if (mattype == "boundary" && matdimension == dimension-1)
        {
            fGmsh.fPZMaterialId[matdimension][matname] = matidcounter;
            InsertDarcyBCMaterial(matidcounter, matdimension, (int)(density+0.5), perm);
            fMaterialIds[matname] = matidcounter;
            matidcounter++;
        }
    }
    ReadNextLine(input, line);
    
    int numstepdefinitions = 0;
    {
        std::istringstream stin(line);
        stin >> numstepdefinitions;
    }
    fTimeSteps.Resize(numstepdefinitions);
    for (int ist=0; ist<numstepdefinitions; ist++) {
        ReadNextLine(input, line);
        std::istringstream stin(line);
        stin >> fTimeSteps[ist].first >> fTimeSteps[ist].second;
    }
    
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> fPostProcessRootname;
    }
    
    int numpostprocess;
    ReadNextLine(input, line);
    {
        std::istringstream stin(line);
        stin >> numpostprocess;
    }
    fPostProcnames.Resize(numpostprocess);
    for (int pp=0; pp<numpostprocess; pp++) {
        ReadNextLine(input, line);
        std::istringstream stin(line);
        stin >> fPostProcnames[pp].first >> fPostProcnames[pp].second;
        fPostProcnames[pp].first = stripquotes(fPostProcnames[pp].first);
        std::cout << "first " << fPostProcnames[pp].first << " second " << fPostProcnames[pp].second << std::endl;
    }

    ReadNextLine(input, line);
    while (line.find("end") == std::string::npos) {
        ReadNextLine(input, line);
        if (!input) {
            std::cout << "No end found in the preamble of the file\n";
            DebugStop();
        }
    }
}



#include "pzvec_extras.h"

TPZTransform<REAL> GetTransform2d(std::set<int64_t> &nodes, TPZGeoMesh *gmesh)
{
    TPZTransform<REAL> result(1,3);
    TPZManVector<int64_t,5> nodevec(nodes.size());
    int i=0;
    for (auto it = nodes.begin() ; it != nodes.end(); it++) {
        nodevec[i++] = *it;
    }
    TPZFNMatrix<25> dist(nodes.size(),nodes.size(),0.);
    for (int i=0; i<nodevec.size(); i++) {
        for (int j=0; j<nodevec.size(); j++) {
            dist(i,j) = 0;
            int64_t in = nodevec[i];
            int64_t jn = nodevec[j];
            for (int c=0; c<3; c++) {
                dist(i,j) += (gmesh->NodeVec()[in].Coord(c)-gmesh->NodeVec()[jn].Coord(c))*
                (gmesh->NodeVec()[in].Coord(c)-gmesh->NodeVec()[jn].Coord(c));
            }
            dist(i,j) = sqrt(dist(i,j));
        }
    }
    REAL maxdist = 0.;
    int jmax = 0;
    int imax = 0;
    for (int i=0; i<nodevec.size(); i++) {
        for (int j=0; j<nodevec.size(); j++) {
            if (dist(i,j) > maxdist) {
                maxdist = dist(i,j);
                imax = i;
                jmax = j;
            }
        }
    }
    TPZManVector<REAL,3> ivec(3),jvec(3),vector(3);
    gmesh->NodeVec()[nodevec[imax]].GetCoordinates(ivec);
    gmesh->NodeVec()[nodevec[jmax]].GetCoordinates(jvec);
    vector = jvec-ivec;
    REAL sum = 0.;
    for (int i=0; i<3; i++) {
        result.Mult()(0,i) = vector[i]*2./(maxdist*maxdist);
        sum += ivec[i]*result.Mult()(0,i);
    }
    result.Sum()(0,0) = -1.-sum;
#ifdef PZDEBUG
    TPZManVector<REAL,2> itr(1),jtr(1);
    result.Apply(ivec, itr);
    result.Apply(jvec, jtr);
    if (abs(itr[0]+1.) > 1.e-6 || abs(jtr[0]-1.) > 1.e-6) {
        DebugStop();
    }
#endif
    return result;
}

TPZTransform<REAL> GetTransform3d(std::set<int64_t> &nodes, TPZGeoMesh *gmesh)
{
    TPZTransform<REAL> result(2,3);
    TPZManVector<REAL,3> minx(3), maxx(3);
    TPZManVector<int64_t,5> nodevec(nodes.size());
    int i=0;
    for (auto it = nodes.begin() ; it != nodes.end(); it++) {
        nodevec[i++] = *it;
    }
    int64_t node0 = nodevec[0];
    gmesh->NodeVec()[node0].GetCoordinates(minx);
    maxx = minx;
    
    for (int i=0; i<nodevec.size(); i++) {
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[nodevec[i]].GetCoordinates(xco);
        for (int ic=0; ic<3; ic++) {
            if(minx[ic] > xco[ic]) minx[ic] = xco[ic];
            if (maxx[ic] < xco[ic]) {
                maxx[ic] = xco[ic];
            }
        }
    }
    if (maxx[0]- minx[0] < 1.e-6) {
        // the surface is perpendicular to the x coordinate
        REAL delx = maxx[1]-minx[1];
        REAL dely = maxx[2]-minx[2];
        result.Mult()(0,1) = 1./(delx/2.);
        result.Mult()(1,2) = 1./(dely/2.);
        TPZManVector<REAL,2> center(2,0.);
        for (i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                center[i] += result.Mult()(i,j)*((minx[j]+maxx[j])*0.5);
            }
        }
        result.Sum()(0,0) = -center[0];
        result.Sum()(1,0) = -center[1];
    }
    else if (maxx[1]- minx[1] < 1.e-6) {
        // the surface is perpendicular to the x coordinate
        REAL delx = maxx[0]-minx[0];
        REAL dely = maxx[2]-minx[2];
        result.Mult()(0,0) = 1./(delx/2.);
        result.Mult()(1,2) = 1./(dely/2.);
        TPZManVector<REAL,2> center(2,0.);
        for (i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                center[i] += result.Mult()(i,j)*(minx[j]+maxx[j])*0.5;
            }
        }
        result.Sum()(0,0) = -center[0];
        result.Sum()(1,0) = -center[1];
    }
    else if (maxx[2]- minx[2] < 1.e-6) {
        // the surface is perpendicular to the x coordinate
        REAL delx = maxx[0]-minx[0];
        REAL dely = maxx[1]-minx[1];
        result.Mult()(0,0) = 1./(delx/2.);
        result.Mult()(1,1) = 1./(dely/2.);
        TPZManVector<REAL,2> center(2,0.);
        for (i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                center[i] += result.Mult()(i,j)*(minx[j]+maxx[j])*0.5;
            }
        }
        result.Sum()(0,0) = -center[0];
        result.Sum()(1,0) = -center[1];
    }
    else
    {
        DebugStop();
    }
        
#ifdef PZDEBUG
    TPZManVector<REAL,2> itr(2),jtr(2);
    result.Apply(minx, itr);
    result.Apply(maxx, jtr);
    for (int i=0; i<2; i++)
    {
        if (abs(itr[i]+1.) > 1.e-6 || abs(jtr[i]-1.) > 1.e-6) {
            DebugStop();
        }
    }
#endif
    return result;
}

void IdentifyCorners2d(std::set<int64_t> &nodeids, TPZGeoMesh *gmesh, TPZTransform<REAL> &tr, TPZVec<int64_t> &nodes)
{
    REAL minco(0.);
    REAL maxco(0.);
    nodes.Resize(2, 0);

    for (auto it=nodeids.begin(); it != nodeids.end(); it++) {
        TPZManVector<REAL,3> coord(3), cotr(1);
        gmesh->NodeVec()[*it].GetCoordinates(coord);
        tr.Apply(coord, cotr);

        if (minco > cotr[0]) {
            minco = cotr[0];
            nodes[0] = *it;
        }
        if (maxco < cotr[0]) {
            maxco = cotr[0];
            nodes[1] = *it;
        }
    }
}

void IdentifyCorners3d(std::set<int64_t> &nodeids, TPZGeoMesh *gmesh, TPZTransform<REAL> &tr, TPZVec<int64_t> &nodes)
{
    int meshdim = gmesh->Dimension();
    if (meshdim != 3) {
        DebugStop();
    }
    nodes.Resize(4, -1);
    nodes.Fill(-1);
    for (auto it=nodeids.begin(); it != nodeids.end(); it++) {
        TPZManVector<REAL,3> coord(3), cotr(meshdim-1);
        gmesh->NodeVec()[*it].GetCoordinates(coord);
        tr.Apply(coord, cotr);
        // lowerleft
        if (fabs(cotr[0]+1.) < 1.e-6 && fabs(cotr[1]+1.) < 1.e-6) {
            nodes[0] = *it;
        }
        // lowerright
        if (fabs(cotr[0]-1.) < 1.e-6 && fabs(cotr[1]+1.) < 1.e-6) {
            nodes[1] = *it;
        }
        // upperright
        if (fabs(cotr[0]-1.) < 1.e-6 && fabs(cotr[1]-1.) < 1.e-6) {
            nodes[2] = *it;
        }
        // upperleft
        if (fabs(cotr[0]+1.) < 1.e-6 && fabs(cotr[1]-1.) < 1.e-6) {
            nodes[3] = *it;
        }
    }
    for (int i=0; i<4; i++) {
        if(nodes[i] == -1)
        {
            DebugStop();
        }
    }
}

int64_t GroupElements(TPZStack<int64_t> &elems, TPZGeoMesh *gmesh)
{
    std::set<int64_t> nodes;
    for (int el=0; el<elems.size(); el++) {
        TPZGeoEl *gel = gmesh->Element(elems[el]);
        int nnodes = gel->NCornerNodes();
        for (int in=0; in<nnodes; in++) {
            nodes.insert(gel->NodeIndex(in));
        }
    }
    int meshdim = gmesh->Dimension();
    TPZTransform<REAL> tr;
    // this method will compute
    if (meshdim == 2) {
        tr = GetTransform2d(nodes, gmesh);
    }
    else if(meshdim == 3)
    {
        tr = GetTransform3d(nodes, gmesh);
    }
    else
    {
        DebugStop();
    }
    
    TPZManVector<int64_t,4> cornernodes;
    if (meshdim == 2) {
        IdentifyCorners2d(nodes, gmesh, tr, cornernodes);
    }
    else if(meshdim == 3)
    {
        IdentifyCorners3d(nodes, gmesh, tr, cornernodes);
    }
    std::map<int64_t,int> nodemap;
    {
        int i = 0;
        for (auto it=nodes.begin(); it != nodes.end(); it++) {
            nodemap[*it] = i++;
        }
    }
    TPZAutoPointer<TPZRefPattern> refpat;
    {
        TPZGeoMesh refpatmesh;
        refpatmesh.NodeVec().Resize(nodemap.size());
        for (auto it=nodemap.begin(); it != nodemap.end(); it++)
        {
            TPZManVector<REAL,3> co(3), cotr(meshdim-1);
            gmesh->NodeVec()[it->first].GetCoordinates(co);
            tr.Apply(co, cotr);
            cotr.Resize(3, 0.);
            refpatmesh.NodeVec()[it->second].Initialize(cotr, refpatmesh);
        }
        TPZManVector<int64_t,4> nodeindices(cornernodes);
        for (int i=0; i<cornernodes.size(); i++) {
            nodeindices[i] = nodemap[cornernodes[i]];
        }
        int64_t index;
        if (meshdim == 2) {
            refpatmesh.CreateGeoElement(EOned, nodeindices, 1, index);
        }
        else if(meshdim == 3)
        {
            refpatmesh.CreateGeoElement(EQuadrilateral, nodeindices, 1, index);
        }
        for (int64_t el=0; el<elems.size(); el++) {
            TPZGeoEl *gel = gmesh->Element(elems[el]);
            int ncorner = gel->NCornerNodes();
            for (int ic=0; ic<ncorner; ic++) {
                nodeindices[ic] = nodemap[gel->NodeIndex(ic)];
            }
            refpatmesh.CreateGeoElement(gel->Type(), nodeindices, 1, index);
            TPZGeoEl *gelsub = refpatmesh.Element(index);
            gelsub->SetFather((int64_t)0);
        }
        refpatmesh.BuildConnectivity();
        refpat = new TPZRefPattern(refpatmesh);
        for (int64_t el=0; el<refpatmesh.NElements(); el++) {
            refpatmesh.Element(el)->SetFather(-1);
        }
    }
    int64_t index;
    int matid = gmesh->Element(elems[0])->MaterialId();
    if (meshdim == 2) {
        gmesh->CreateGeoElement(EOned, cornernodes, matid, index);
    }
    else
    {
        gmesh->CreateGeoElement(EQuadrilateral, cornernodes, matid, index);
    }
    refpat->GenerateSideRefPatterns();
    
    gmesh->Element(index)->SetRefPattern(refpat);
    for (int el = 0; el<elems.size(); el++) {
        gmesh->Element(index)->SetSubElement(el, gmesh->Element(elems[el]));
        gmesh->Element(elems[el])->SetFather(gmesh->Element(index));
    }
    return index;
}

void CreateRefPatterns(TPZGeoMesh *gmesh, TPZVec<int64_t> &elemententity, std::set<int> &matids)
{
    int64_t nel = gmesh->NElements();
    int meshdim = gmesh->Dimension();
    int64_t el =0;
    while(el<nel)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        // find an element of dimension-1
        if(!gel || (gel->Dimension() != meshdim-1))
        {
            el++;
            continue;
        }
        int matid = gel->MaterialId();
        if (matids.find(matid) == matids.end()) {
            el++;
            continue;
        }
        TPZStack<int64_t> elstack;
        elstack.Push(el);
        int64_t ident = elemententity[el];
        el++;
        while(el < nel && elemententity[el] == ident)
        {
            elstack.Push(el);
            el++;
        }
        if (elstack.size() > 1) {
            int64_t index = GroupElements(elstack, gmesh);
            if (elemententity.size() <= index) {
                elemententity.Resize(index+1, -1);
            }
            elemententity[index] = ident;
        }
    }
    gmesh->BuildConnectivity();
}

/// verify is each macro element links only two subdomains
void VerifySubdomainConsistency(TPZGeoMesh *gmesh, TPZVec<int64_t> &elemententity, std::set<int> &matids)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if (gel->Dimension() != gmesh->Dimension()-1 || matids.find(matid) == matids.end()) {
            continue;
        }
        if (!gel->HasSubElement()) {
            continue;
        }
        std::set<int64_t> connected;
        int nsubel = gel->NSubElements();
        for (int isub = 0; isub<nsubel; isub++) {
            TPZGeoEl *subel = gel->SubElement(isub);
            TPZGeoElSide subelside(subel,subel->NSides()-1);
            TPZGeoElSide neighbour = subelside.Neighbour();
            while (neighbour != subelside) {
                if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
                    connected.insert(elemententity[neighbour.Element()->Index()]);
                }
                neighbour = neighbour.Neighbour();
            }
        }
        if (connected.size() != 2) {
            std::cout << "Coordinates of the element\n";
            TPZManVector<REAL,3> co0(3), co1(3);
            gel->Node(0).GetCoordinates(co0);
            gel->Node(1).GetCoordinates(co1);
            std::cout << co0 << " " << co1 << std::endl;
            DebugStop();
        }
    }
}

/// Build skeleton data structure
void BuildSkeleton(TPZGeoMesh *gmesh, TPZVec<int64_t> &elemententity, std::set<int> &matids, std::map<int64_t,std::pair<int64_t,int64_t>> &skeleton)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if (gel->Dimension() != gmesh->Dimension()-1 || matids.find(matid) == matids.end()) {
            continue;
        }
        if (gel->Father()) {
            continue;
        }
        std::set<int64_t> connected;
        TPZStack<TPZGeoElSide> elstack;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        elstack.Push(gelside);
        while(elstack.size())
        {
            TPZGeoElSide search = elstack.Pop();
            TPZGeoElSide neighbour = search.Neighbour();
            while (neighbour != search)
            {
                if (neighbour.Element()->Dimension() == gmesh->Dimension())
                {
                    connected.insert(elemententity[neighbour.Element()->Index()]);
                }
                neighbour = neighbour.Neighbour();
            }
            if (search.Element()->HasSubElement())
            {
                int nsubel = search.Element()->NSubElements();
                for (int isub = 0; isub<nsubel; isub++)
                {
                    TPZGeoEl *subel = search.Element()->SubElement(isub);
                    TPZGeoElSide subelside(subel,subel->NSides()-1);
                    elstack.Push(subelside);
                }
            }
        }
#ifdef PZDEBUG
        if (connected.size() > 2)
        {
            std::cout << "Coordinates of the element\n";
            TPZManVector<REAL,3> co0(3), co1(3);
            gel->Node(0).GetCoordinates(co0);
            gel->Node(1).GetCoordinates(co1);
            std::cout << co0 << " " << co1 << std::endl;
            DebugStop();
        }
#endif
        if(connected.size() == 2)
        {
            auto it = connected.begin();
            int64_t domain1 = *it;
            it++;
            int64_t domain2 = *it;
            skeleton[el] = std::pair<int64_t,int64_t>(domain1,domain2);
        } else if (connected.size() == 1)
        {
            auto it = connected.begin();
            int64_t domain1 = *it;
            int64_t domain2 = el;
            skeleton[el] = std::pair<int64_t, int64_t>(domain1,domain2);
        }
        else
        {
            DebugStop();
        }
    }
}

/// Adjust the entity index of the fracture elemenents
void AdjustEntityOfFractures(TPZGeoMesh *gmesh,TPZVec<int64_t> &EntityIndex, std::set<int> fracmatid)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (fracmatid.find(gel->MaterialId()) != fracmatid.end()) {
            std::set<int64_t> entities;
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == gmesh->Dimension()) {
                    entities.insert(EntityIndex[neighbour.Element()->Index()]);
                }
                neighbour = neighbour.Neighbour();
            }
            if (entities.size() != 1) {
                TPZManVector<REAL,3> co0(3),co1(3);
                gel->Node(0).GetCoordinates(co0);
                gel->Node(1).GetCoordinates(co1);
                std::cout << "Coordinates of the element " << co0 << " " << co1 << std::endl;
                DebugStop();
            }
            EntityIndex[gel->Index()] = *entities.begin();
        }
    }
}


/// creates and inserts a Darcy or ParabolicDarcy object in the mesh
void TPZFracSimulation::InsertDarcyMaterial(int matid, REAL permeability, REAL rho)
{
#ifdef PZDEBUG
    if (fMHM->fMaterialIds.find(matid) != fMHM->fMaterialIds.end()) {
        DebugStop();
    }
#endif
    int dimension = fMHM->GMesh()->Dimension();
    TPZMixedPoisson * mat;
    if (fSimulationType == 0)
    {
        mat = new TPZMixedPoisson(matid,dimension);
    }
    else
    {
        TPZMixedPoissonParabolic *matp = new TPZMixedPoissonParabolic(matid,dimension);
        matp->SetDeltaT(1000.);
        mat = matp;
    }
    mat->SetSymmetric();
    mat->SetPermeability(permeability);
    
    fMHM->CMesh()->InsertMaterialObject(mat);
    
//    TPZVecL2 *vecmat = new TPZVecL2(matid);
//    fMHM->FluxMesh()->InsertMaterialObject(vecmat);
//    TPZMatLaplacian *presmat = new TPZMatLaplacian(matid);
//    presmat->SetDimension(dimension);
//    fMHM->PressureMesh()->InsertMaterialObject(presmat);
    
    fMHM->fMaterialIds.insert(matid);
    
}

/// creates and inserts the boundary condition objects
void TPZFracSimulation::InsertDarcyBCMaterial(int matid, int dimension, int bctype, REAL val)
{
#ifdef PZDEBUG
    if (fMHM->fMaterialBCIds.find(matid) != fMHM->fMaterialBCIds.end()) {
        DebugStop();
    }
#endif
    int rootmat = *fMHM->fMaterialIds.begin();
    TPZMaterial *mat = fMHM->CMesh()->FindMaterial(rootmat);
    if (!mat) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,val);
    TPZBndCond * bcmat = new TPZBndCond(mat,matid,bctype,val1,val2);
    //    TPZMixedPoissonParabolic *mat = new TPZMixedPoissonParabolic(matid,dimension);
    //    mat->SetDeltaT(1000.);
    
    fMHM->CMesh()->InsertMaterialObject(bcmat);
    
    
    fMHM->fMaterialBCIds.insert(matid);
    
}


void TPZFracSimulation::ReadFractures(std::ifstream &input)
{
    std::string line;
    int dimension = fMHM->GMesh()->Dimension();
    int id = 0;
    while(input)
    {
        if(std::getline(input, line))
        {
            if (line[0] == '#') {
                continue;
            }
            if(line.find("CORNER") != std::string::npos && dimension == 2)
            {
                TPZGeoNode first, second;
                TPZManVector<REAL,3> x0(3), x1(3);
                std::istringstream sin(line);
                std::string corner;
                sin >> corner >> x0[0] >> x0[1] >> x0[2] >> x1[0] >> x1[1] >> x1[2];
                first.SetCoord(x0);
                second.SetCoord(x1);
                int64_t index0 = fFracSet.InsertNode(first);
                int64_t index1 = fFracSet.InsertNode(second);
                first = fFracSet.fNodeVec[index0];
                second = fFracSet.fNodeVec[index1];
                uint64_t keyfirst = fFracSet.GetLoc(first);
                uint64_t keysecond = fFracSet.GetLoc(second);
                if (fFracSet.fPointMap.find(keyfirst) == fFracSet.fPointMap.end()) {
                    DebugStop();
                }
                if (fFracSet.GetLoc(first) != keyfirst) {
                    DebugStop();
                }
                if (fFracSet.fPointMap.find(keysecond) == fFracSet.fPointMap.end()) {
                    DebugStop();
                }
                if (fFracSet.GetLoc(second) != keysecond) {
                    DebugStop();
                }
                int64_t lastfrac = fFracSet.fFractureVec.NElements()-1;
                if (lastfrac < 0) {
                    DebugStop();
                }

                if (first.Coord(0) < second.Coord(0)) {
                    int matid = fFracSet.matid_internal_frac;
                    TPZFracture frac(id,matid,index0,index1);
                    frac.fPhysicalName = fFracSet.fFractureVec[lastfrac].fPhysicalName;
                    fFracSet.fFractureVec[lastfrac] = frac;
                }
                else
                {
                    int matid = fFracSet.matid_internal_frac;
                    TPZFracture frac(id,matid,index1,index0);
                    frac.fPhysicalName = fFracSet.fFractureVec[lastfrac].fPhysicalName;
                    fFracSet.fFractureVec[lastfrac] = frac;
                }
            }
            {
                uint64_t pos = line.find("PERM");
                if(pos != std::string::npos)
                {
                    int64_t lastfrac = fFracSet.fFractureVec.NElements()-1;
                    if(lastfrac < 0) DebugStop();
                    std::string sub = line.substr(pos+4,std::string::npos);
                    std::istringstream sin(sub);
                    sin >> fFracSet.fFractureVec[lastfrac].fFracPerm;
                }
            }
            {
                uint64_t pos = line.find("THICK");
                if(pos != std::string::npos)
                {
                    int64_t lastfrac = fFracSet.fFractureVec.NElements()-1;
                    if(lastfrac < 0) DebugStop();
                    std::string sub = line.substr(pos+5,std::string::npos);
                    std::istringstream sin(sub);
                    double thickness;
                    sin >> thickness;
                    fFracSet.fFractureVec[lastfrac].fThickness = thickness;
                }
            }
            {
                uint64_t pos = line.find("NAME");
                if(pos != std::string::npos)
                {
                    int64_t index = fFracSet.fFractureVec.AllocateNewElement();
                    std::string sub = line.substr(pos+4,std::string::npos);
                    std::istringstream sin(sub);
                    std::string fracname;
                    sin >> fracname;
                    fracname.erase(0,1);
                    fracname.erase(fracname.end()-1,fracname.end());
                    fFracSet.fFractureVec[index].fPhysicalName = fracname;
                }
            }
        }
        else
        {
            break;
        }
    }
    int64_t nfrac = fFracSet.fFractureVec.NElements();
    int matid = *(fMHM->fMaterialBCIds.rbegin())+1;
    matid = (matid+10)-matid%10;
    for (int ifr=0; ifr<nfrac; ifr++) {
        int meshdim = fMHM->GMesh()->Dimension();
        fMHM->InsertFractureFlowMaterial(matid);
        // Material medio poroso
        TPZMixedPoisson * mat = new TPZMixedPoisson(matid,meshdim-1);
        TPZFracture frac = fFracSet.fFractureVec[ifr];
        REAL perm = frac.fFracPerm * frac.fThickness;
        std::string matname = fFracSet.fFractureVec[ifr].fPhysicalName;
        fGmsh.fPZMaterialId[meshdim-1][matname] = matid;
        fMaterialIds[matname] = matid;
        mat->SetSymmetric();
        mat->SetPermeability(perm);
        TPZFNMatrix<9,REAL> K(3,3,0.),KInv(3,3,0.);
        K(0,0) = perm;
        K(1,1) = perm;
        K(2,2) = perm;
        KInv(0,0) = 1./K(0,0);
        KInv(1,1) = 1./K(1,1);
        KInv(2,2) = 1./K(2,2);
        mat->SetPermeabilityTensor(K, KInv);
        fMHM->CMesh()->InsertMaterialObject(mat);
        TPZVecL2 *vecl21 = new TPZVecL2(matid);
        fMHM->FluxMesh()->InsertMaterialObject(vecl21);
//        TPZMat1dLin *mat1d = new TPZMat1dLin(matid);
//        fMHM->PressureMesh()->InsertMaterialObject(mat1d);
        matname = matname + "_MHM";
        matid++;
        mat = new TPZMixedPoisson(*mat);
        mat->SetId(matid);
        fGmsh.fPZMaterialId[meshdim-1][matname] = matid;
        fMaterialIds[matname] = matid;
        fMHM->CMesh()->InsertMaterialObject(mat);
        TPZVecL2 *vecl2 = new TPZVecL2(matid);
        fMHM->FluxMesh()->InsertMaterialObject(vecl2);
//        TPZMat1dLin *mat1d2 = new TPZMat1dLin(matid);
//        fMHM->PressureMesh()->InsertMaterialObject(mat1d2);

        fMHM->fSkeletonWithFlowMatId.insert(matid);
        matid++;
    }
}

/// adjust the geometric element read from gmesh
void TPZFracSimulation::AdjustGeometricMesh(const std::string &rootname)
{
    {
        std::string filename = rootname + "_gmesh.vtk";
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh().operator->(), out);
    }
    
    TPZGeoMesh *gmesh = fMHM->GMesh().operator->();
    std::set<int> matids = fMHM->fSkeletonWithFlowMatId;
    matids.insert(fMHM->fSkeletonMatId);
    CreateRefPatterns(gmesh, fGmsh.fEntityIndex, matids);
    VerifySubdomainConsistency(gmesh, fGmsh.fEntityIndex, matids);
    matids = fMHM->fFractureFlowDim1MatId;
    AdjustEntityOfFractures(gmesh,fGmsh.fEntityIndex,matids);
    std::map<int64_t,std::pair<int64_t,int64_t>> skeletonstruct;
    matids = fMHM->fSkeletonWithFlowMatId;
    matids.insert(fMHM->fSkeletonMatId);
    std::copy(fMHM->fMaterialBCIds.begin(),fMHM->fMaterialBCIds.end(), std::inserter(matids, matids.begin()));
    BuildSkeleton(gmesh, fGmsh.fEntityIndex, matids, skeletonstruct);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        sout << "Entity index \n" << fGmsh.fEntityIndex << endl;
        for (auto it = skeletonstruct.begin(); it != skeletonstruct.end(); it++) {
            TPZGeoEl *gel = fMHM->GMesh()->Element(it->first);
            sout << "skeleton geo " << it->first << " type " << gel->TypeName() << " mat id " << gel->MaterialId() << " linked to subdomains " << it->second.first << " " << it->second.second;
            if (it->first == it->second.second) {
                sout << " boundary";
            }
            sout << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef PZDEBUG
    if(0)
    {
        int numquad = 0;
        std::set<int> matids;
        for (int64_t el=0; el< gmesh->NElements(); el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel && gel->Type() == EQuadrilateral)
            {
//                std::cout << gel->SideArea(8) << std::endl;
                matids.insert(gel->MaterialId());
                numquad++;
            }
        }
        std::cout << "number of quadrilaterals = " << numquad << std::endl;
    }
#endif
#ifdef PZDEBUG
    {
        std::string filename = rootname + ".gmesh_mhm.vtk";
        std::ofstream out(filename);
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh(), out);
    }
    {
        std::string filename = rootname + ".gmesh_eltype.vtk";
        std::ofstream out(filename);
        TPZVec<int> eltype(fMHM->GMesh()->NElements(),-1);
        for (int64_t el=0; el<fMHM->GMesh()->NElements(); el++) {
            TPZGeoEl *gel = fMHM->GMesh()->Element(el);
            if (gel) {
                eltype[el] = gel->Type();
            }
        }
        TPZVTKGeoMesh::PrintGMeshVTK(fMHM->GMesh().operator->(), out, eltype);
    }
#endif
    fMHM->DefinePartition(fGmsh.fEntityIndex, skeletonstruct);
}

/// adjust the boundary condition type of the elements
void TPZFracSimulation::SetStandardBoundaryElements()
{
    TPZGeoMesh *gmesh = fMHM->GMesh().operator->();
    int dimension = gmesh->Dimension();
    int matidbc = fGmsh.fPZMaterialId[dimension-1]["BC"];
    int64_t nel = gmesh->NElements();
    for (int64_t el =0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel && gel->Dimension() == dimension-1 && gel->MaterialId() == matidbc) {
            TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
            int nnodes = gel->NCornerNodes();
            for (int in = 0; in<nnodes; in++) {
                TPZManVector<REAL, 3> co(3);
                gel->Node(in).GetCoordinates(co);
                if (in == 0) {
                    minx = co;
                    maxx = co;
                }
                for (int ic=0 ; ic<3; ic++) {
                    if (minx[ic] < co[ic]) {
                        minx[ic] = co[ic];
                    }
                    if (maxx[ic] > co[ic]) {
                        maxx[ic] = co[ic];
                    }
                }
            }
            if (fabs(minx[0]-maxx[0]) < 1.e-6 && fabs(minx[0]-fFracSet.fLowLeft[0]) < 1.e-6) {
                int matid = fGmsh.fPZMaterialId[dimension-1]["BCIN"];
                gel->SetMaterialId(matid);
            }
            else if (fabs(minx[0]-maxx[0]) < 1.e-6 && fabs(maxx[0]-fFracSet.fTopRight[0]) < 1.e-6) {
                int matid = fGmsh.fPZMaterialId[dimension-1]["BCOUT"];
                gel->SetMaterialId(matid);
            }
            // top or bottom
            else {
                int matid = fGmsh.fPZMaterialId[dimension-1]["BCNOFLOW"];
                gel->SetMaterialId(matid);
            }
        }
    }
}


