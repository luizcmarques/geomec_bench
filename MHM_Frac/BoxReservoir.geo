// Gmsh project created on Sun Jun 24 10:44:31 2018
SetFactory("OpenCASCADE");
//+
Block(1) = {0, 0, 0, 50, 50, 10};

Physical Volume("reservoir") = {1};
Physical Surface("top") = {5,6};
Physical Surface("inflow") = {1};
Physical Surface("lateral") = {3,4};
Physical Surface("outflow") = {2};



Transfinite Line {2,4,6,8,9,10,11,12} = 11 Using Progression 1;
Transfinite Line {1,3,5,7} = 6 Using Progression 1;


Transfinite Surface {1,2,3,4,5,6};


Recombine Surface {1,2,3,4,5,6};


Transfinite Volume {1};

