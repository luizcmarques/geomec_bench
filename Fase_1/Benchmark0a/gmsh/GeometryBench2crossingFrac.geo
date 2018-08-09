
SetFactory("OpenCASCADE");

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Parameters definitions

a = 5; 
b = 5;
h = 10;
L = 200;
Lf = 200;

n_bc = 2;
nx = 10;
ny = 4;
n_f[] = {20,4};
pr = 1.0;

  // Points

  //Omega
  p1 = newp; Point(p1) = {0, 0, 0};
  p2 = newp; Point(p2) = {L, 0, 0};
  p3 = newp; Point(p3) = {L, h, 0};
  p4 = newp; Point(p4) = {0, h, 0};

  //Points for fracture 1
  p5 = newp; Point(p5) = {0, a, 0};
  p6 = newp; Point(p6) = {150, b, 0};  

  //Points for fracture 2
  p7 = newp; Point(p7) = {145, 7, 0};
  p8 = newp; Point(p8) = {155, 3, 0}; 

  // Boundaries
  l1 = newl; Line(l1) = {p1, p2};
  l2 = newl; Line(l2) = {p2, p3};
  l3 = newl; Line(l3) = {p3, p4};
  l4 = newl; Line(l4) = {p4, p5};  
  l5 = newl; Line(l5) = {p5, p1};  

  // Frature 1
  f1 = newl; Line(f1) = {p5, p6};

  // Frature 2
  f2 = newl; Line(f2) = {p7, p8};

  // f1 -> f2
  binary_fracture_intersection[] = BooleanFragments{ Line{f1}; }{ Line{f2}; }; 
  Printf("Intersection created the lines ", binary_fracture_intersection[]);

  n_lines = #binary_fracture_intersection[];
  If(n_lines==2)
    Printf("No intersection with provided lines ",{f1,f2});
    f1[] = {binary_fracture_intersection[0]};
    f2[] = {binary_fracture_intersection[1]};
  Else
    f1[] = {binary_fracture_intersection[0],binary_fracture_intersection[1]};
    f2[] = {binary_fracture_intersection[2],binary_fracture_intersection[3]};
  EndIf  

  // Refinements

  Transfinite Line{l1,l3} = nx Using Progression pr;
  Transfinite Line{l2,l4} = ny Using Progression pr;


  Transfinite Line{f1} = n_f[0] Using Progression pr;
  Transfinite Line{f2} = n_f[1] Using Progression pr;
  //Transfinite Line{f3} = n_f[2] Using Progression pr;

// Definição da superfície 
  ll1 = newll; Line Loop(ll1) = {l1, l2, l3, l4, l5};
  s1 = news; Plane Surface(s1) = {ll1};

// Intersection point



// Associating fractures and fracture boundaries to Omega
  fracture_set[] = {f2,f1};
  fracture_point_set[] = {p5,p6,p7,p8};
  Line{fracture_set[]} In Surface{s1};
  Point{fracture_point_set[]} In Surface{s1};

// Computing the intersection and clustering intersected fractures
  disjunctive_p = p6;
  fracture_cluster[] = fracture_set[];

  Point{disjunctive_p} In Line{fracture_cluster[]};

 // Transfinite Surface {s1};

  If(IsquadQ)
    Recombine Surface {s1};
  EndIf


  Physical Surface("Omega") = {s1};
  Physical Line("bottom") = {l1};
  Physical Line("right") = {l2};
  Physical Line("top") = {l3};
  Physical Line("left") = {l4};
  
  Physical Line("frac") = {f1};
  Physical Point("PointLeft") = {p5};
  Physical Point("PointRight") = {p6};

  Physical Line("frac1") = {f2};
  Physical Point("PointLeft1") = {p7};
  //Physical Point("PointRight1") = {p6};

  //Physical Line("frac2") = {f3};
  //Physical Point("PointLeft2") = {p6};
  Physical Point("PointRight2") = {p8};
  
  Coherence Mesh;


