

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

a = 5; 
b = 5;
h = 10;
L = 200;
Lf = 200;

n_bc = 2;
nx = 12;
ny = 2;
pr = 1.0;

// Coordenadas dos pontos

  //Domínio Omega
  p1 = newp; Point(p1) = {0, 0, 0};
  p2 = newp; Point(p2) = {L, 0, 0};
  p3 = newp; Point(p3) = {L, h, 0};
  p4 = newp; Point(p4) = {0, h, 0};

  //Fratura
  p5 = newp; Point(p5) = {50, a, 0};
  p6 = newp; Point(p6) = {150, b, 0};  

// Fronteiras

  //Domínio Omega  
  l1 = newl; Line(l1) = {p1, p2};
  l2 = newl; Line(l2) = {p2, p3};
  l3 = newl; Line(l3) = {p3, p4};
  l4 = newl; Line(l4) = {p4, p1};  

  //Fratura
  f1 = newl; Line(f1) = {p5, p6};

  Transfinite Line{l1,l3} = nx Using Progression pr;
  Transfinite Line{l2,l4} = ny Using Progression pr;
  //Transfinite Line{f1} = nx Using Progression pr;


// Definição da superfície 
  Line Loop(1) = {l1, l2, l3, l4};
  Plane Surface(1) = {1};
  Line{f1} In Surface{1};

 // Point{p5,p6} In Surface{s1};

 // Transfinite Surface {s1};

  If(IsquadQ)
    Recombine Surface {1};
  EndIf


  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {l1};
  Physical Line("right") = {l2};
  Physical Line("top") = {l3};
  Physical Line("left") = {l4};

  Physical Line("frac") = {f1};
  
  Coherence Mesh;


