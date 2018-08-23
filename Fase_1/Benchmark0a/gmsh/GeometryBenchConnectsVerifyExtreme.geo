

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

a = 5; 
b = 5;
h = 10;
L = 10;
Lf = 10;

n_bc = 2;
nx = 2;
ny = 2;
pr = 1.0;

// Coordenadas dos pontos

  //Domínio Omega
  p1 = newp; Point(p1) = {0, 0, 0};
  p2 = newp; Point(p2) = {L, 0, 0};
  p3 = newp; Point(p3) = {L, h, 0};
  p4 = newp; Point(p4) = {0, h, 0};

  //Fratura
  p5 = newp; Point(p5) = {0, a, 0};
  p6 = newp; Point(p6) = {10, 6.8, 0};  

// Fronteiras

  //Domínio Omega  
  l1 = newl; Line(l1) = {p1, p2};
  l2 = newl; Line(l2) = {p2, p6};
  l3 = newl; Line(l3) = {p6, p3};
  l4 = newl; Line(l4) = {p3, p4};
  l5 = newl; Line(l5) = {p4, p5};  
  l6 = newl; Line(l6) = {p5, p1};  

  //Fratura
  f1 = newl; Line(f1) = {p5, p6};

  Transfinite Line{l1,l4} = nx Using Progression pr;
  Transfinite Line{l2,l3,l5,l6} = ny Using Progression pr;
  Transfinite Line{f1} = 2 Using Progression pr;


// Definição da superfície 
  ll1 = newll; Line Loop(ll1) = {l1, l2, l3, l4, l5, l6};
  s1 = news; Plane Surface(s1) = {ll1};
  Line{f1} In Surface{s1};
  Point{p5,p6} In Surface{s1};

 // Transfinite Surface {s1};

  If(IsquadQ)
    Recombine Surface {s1};
  EndIf


  Physical Surface("Omega") = {s1};
  Physical Line("bottom") = {l1};
  Physical Line("right") = {l2,l3};
  Physical Line("top") = {l4};
  Physical Line("left") = {l5,l6};
  Physical Line("frac") = {f1};
  Physical Point("PointLeft") = {p5};
  Physical Point("PointRight") = {p6};
  
  Coherence Mesh;


