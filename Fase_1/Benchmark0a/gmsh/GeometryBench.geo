

IsquadQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

a = 5; 
b = 5;
h = 10;
L = 200;
Lf = 200;

n_bc = 2;
ny = 8;
nx = 2;
pr = 1;

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {L, 0, 0, 1e+22};
  Point(3) = {L, h, 0, 1e+22};
  Point(4) = {0, h, 0, 1e+22};

  //Fratura
  Point(5) = {0, a, 0, 1e+22};
  Point(6) = {L, b, 0, 1e+22};  


// Fronteiras

  //Domínio Omega  
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};  

  //Fratura
  Line(5) = {5, 6};

  Transfinite Line{1,3} = ny Using Progression pr;
  Transfinite Line{2,4} = nx Using Progression pr;
  Transfinite Line{5} = ny Using Progression pr;


// Definição da superfície 

  Line Loop(5) = {1, 2, 3, 4};
 
  Plane Surface(6) = {5};

  If(IsquadQ)

  Recombine Surface {6};

  EndIf

 // Physical Volume("internal") = {1};
 // Extrude {0, 0, 10} {
 //  //Surface{110};
 //  Layers{1};
 //  Recombine;
 // }

  Physical Surface("Omega") = {6};
  Physical Line("bottom") = {1};
  Physical Line("top") = {3};
  Physical Line("right") = {2};
  Physical Line("left") = {4};
  Physical Line("frac") = {7};
  
  //Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




