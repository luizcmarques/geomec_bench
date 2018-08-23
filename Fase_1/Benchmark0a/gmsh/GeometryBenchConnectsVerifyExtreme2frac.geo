
SetFactory("OpenCASCADE");

Include "macros/IntersectLines.geo";

IsquadQ = 1;
 
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
n_f[] = {2,2};
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

  //Fratura 2
  p7 = newp; Point(p7) = {5, 0, 0};
  p8 = newp; Point(p8) = {7, 10, 0};  

// Fronteiras

  //Domínio Omega  
  l1 = newl; Line(l1) = {p1, p7};
  l2 = newl; Line(l2) = {p7, p2};
  
  l3 = newl; Line(l3) = {p2, p6};
  l4 = newl; Line(l4) = {p6, p3};
  
  l5 = newl; Line(l5) = {p3, p8};
  l6 = newl; Line(l6) = {p8, p4};
  
  l7 = newl; Line(l7) = {p4, p5};  
  l8 = newl; Line(l8) = {p5, p1};  

  //Fratura
  f1 = newl; Line(f1) = {p5, p6};
  f2 = newl; Line(f2) = {p7, p8};

     bunch_of_entities[] = BooleanFragments{ Line{f1}; }{ Line{f2}; };
    //Printf("intersection created lines ", bunch_of_entities[]);

    point_cluster[] = {};

    n_lines = #bunch_of_entities[];

    If(n_lines == 3)
      f1[] = {bunch_of_entities[0]};
      f2[] = {bunch_of_entities[1],bunch_of_entities[2]};
      
      points[] = Boundary{Line{f1[0]};};
      f1points[] = {points[0],points[1]};

      points_1[] = Boundary{Line{f2[0]};};
      points_2[] = Boundary{Line{f2[1]};};
      f2points[] = {points_1[0],points_2[1]};
    EndIf

    If(n_lines == 4)
      Printf("intersection created lines ", bunch_of_entities[]);
      f1[] = {bunch_of_entities[0],bunch_of_entities[1]};
      f2[] = {bunch_of_entities[2],bunch_of_entities[3]};
      
      points_1[] = Boundary{Line{f1[0]};};
      points_2[] = Boundary{Line{f1[1]};};
      f1points[] = {points_1[0],points_2[1]};

      points_1[] = Boundary{Line{f2[0]};};
      points_2[] = Boundary{Line{f2[1]};};
      f2points[] = {points_1[0],points_2[1]};
    EndIf   

    point_cluster[] = {f1points[],f2points[]};

  Coherence;

  Transfinite Line{l1,l2,l5,l6} = nx Using Progression pr;
  Transfinite Line{l3,l4,l7,l8} = ny Using Progression pr;
//  Transfinite Line{f1} = 2 Using Progression pr;
//  Transfinite Line{f2} = 2 Using Progression pr;
  Transfinite Line{f1[]} = n_f[0] Using Progression pr;
  Transfinite Line{f2[]} = n_f[1] Using Progression pr;

// Definição da superfície 
  ll1 = newll; Line Loop(ll1) = {l1, l2, l3, l4, l5, l6, l7, l8};
  s1 = news; Plane Surface(s1) = {ll1};

// Associating fractures and fracture boundaries to Omega
  fracture_set[] = {f1[],f2[]};
  fracture_point_set[] = {point_cluster[]};
  Line{fracture_set[]} In Surface{s1};

 // Transfinite Surface {s1};

  If(IsquadQ)
    Recombine Surface {s1};
  EndIf


  Physical Surface("Omega") = {s1};
  Physical Line("bottom") = {l1,l2};
  Physical Line("right") = {l3,l4};
  Physical Line("top") = {l5,l6};
  Physical Line("left") = {l7,l8};
  Physical Line("frac") = {f1[]};
  Physical Line("frac2") = {f2[]};
  
  Coherence Mesh;


