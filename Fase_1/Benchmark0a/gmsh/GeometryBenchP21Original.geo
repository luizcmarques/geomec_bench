

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
nx = 18;
ny = 8;
pr = 1;

// Coordenadas dos pontos

  //Domínio Omega
  p1 = newp; Point(p1) = {0, 0, 0, 1e+22};
  p2 = newp; Point(p2) = {L, 0, 0, 1e+22};
  p3 = newp; Point(p3) = {L, h, 0, 1e+22};
  p4 = newp; Point(p4) = {0, h, 0, 1e+22};

  //Fraturas

  p5 = newp; Point(p5) = {59.56005859375, 0, 0, 1e+22}; 
  p6 = newp; Point(p6) = {57.26318359375, 10, 0, 1e+22};

  p7 = newp; Point(p7) = {76.92333984375, 5.39453125, 0, 1e+22};
  p8 = newp; Point(p8) = {76.5068359375, 10, 0, 1e+22};

  p9 = newp; Point(p9) = {89.01513671875, 0, 0, 1e+22};
  p10 = newp; Point(p10) = {87.35107421875, 10, 0, 1e+22};

  p11 = newp; Point(p11) = {97.42724609375,  0, 0, 1e+22};
  p12 = newp; Point(p12) = {90.37890625, 10, 0, 1e+22};

  p13 = newp; Point(p13) = {102.81201171875, 10, 0, 1e+22};
  p14 = newp; Point(p14) = {102.45849609375, 2.873046875, 0, 1e+22};

  p15 = newp; Point(p15) = {127.0439453125,  0, 0, 1e+22};
  p16 = newp; Point(p16) = {121.64697265625, 10, 0, 1e+22};

  p17 = newp; Point(p17) = {191.47509765625, 0, 0, 1e+22};
  p18 = newp; Point(p18) = {188.8212890625,  10, 0, 1e+22};

  p19 = newp; Point(p19) = {157.5068359375,  1.970703125, 0, 1e+22};
  p20 = newp; Point(p20) = {157.3876953125,  0, 0, 1e+22};

  p21 = newp; Point(p21) = {177.08203125,  0, 0, 1e+22};
  p22 = newp; Point(p22) = {173.21044921875, 7.1572265625, 0, 1e+22};

  p23 = newp; Point(p23) = {184.17529296875, 0, 0, 1e+22};
  p24 = newp; Point(p24) = {182.79931640625, 10, 0, 1e+22};

  p25 = newp; Point(p25) = {131.27392578125, 0, 0, 1e+22};
  p26 = newp; Point(p26) = {125.04931640625, 10, 0, 1e+22};

  p27 = newp; Point(p27) = {192.71142578125, 0, 0, 1e+22};
  p28 = newp; Point(p28) = {192.60400390625, 1.470703125, 0, 1e+22};

//  Point(29) = {87.95263671875,  6.6376953125, 0, 1e+22};
//  Point(30) = {87.3271484375, 10, 0, 1e+22};

  p29 = newp; Point(p29) = {193.7646484375,  3.07763671875, 0, 1e+22};
  p30 = newp; Point(p30) = {192.798828125, 10, 0, 1e+22};

  p31 = newp; Point(p31) = {195.78076171875, 0, 0, 1e+22};
  p32 = newp; Point(p32) = {195.13818359375, 1.36083984375, 0, 1e+22};

// Fronteiras

  //Domínio Omega  
  l1 = newl; Line(l1) = {p1,p5};
  l2 = newl; Line(l2) = {p5,p9};
  l3 = newl; Line(l3) = {p9,p11};
  l4 = newl; Line(l4) = {p11,p15};
  l5 = newl; Line(l5) = {p15,p25};  
  l6 = newl; Line(l6) = {p25,p20};
  l7 = newl; Line(l7) = {p20,p21};
  l8 = newl; Line(l8) = {p21,p23};
  l9 = newl; Line(l9) = {p23,p17};
  l10 = newl; Line(l10) = {p17,p27};
  l11 = newl; Line(l11) = {p27,p31};
  l12 = newl; Line(l12) = {p31,p2};
  l13 = newl; Line(l13) = {p2,p3};
  l14 = newl; Line(l14) = {p3,p30};
  l15 = newl; Line(l15) = {p30,p18};
  l16 = newl; Line(l16) = {p18,p24};
  l17 = newl; Line(l17) = {p24,p26};
  l18 = newl; Line(l18) = {p26,p16};
  l19 = newl; Line(l19) = {p16,p13};
  l20 = newl; Line(l20) = {p13,p12};
  l21 = newl; Line(l21) = {p12,p10};
  l22 = newl; Line(l22) = {p10,p8};
  l23 = newl; Line(l23) = {p8,p6};
  l24 = newl; Line(l24) = {p6,p4};
  l25 = newl; Line(l25) = {p4,p1};

  //Fratura

  fshift = 25;

  f1 = newl; Line(f1) = {p5,p6};
  f2 = newl; Line(f2) = {p7,p8};
  f3 = newl; Line(f3) = {p9,p10};
  f4 = newl; Line(f4) = {p11,p12};
  f5 = newl; Line(f5) = {p13,p14};
  f6 = newl; Line(f6) = {p15,p16};
  f7 = newl; Line(f7) = {p17,p18};
  f8 = newl; Line(f8) = {p19,p20};
  f9 = newl; Line(f9) = {p21,p22};
  f10 = newl; Line(f10) = {p23,p24};
  f11 = newl; Line(f11) = {p25,p26};
  f12 = newl; Line(f12) = {p27,p28};
//  Line() = {29, 30};
  f13 = newl; Line(f13) = {p29,p30};
  f14 = newl; Line(f14) = {p31,p32};

  Transfinite Line{l1,l2,l4,l6,l7,l8,l13,l14,l15,l16,l17,l19,l20,l22,l24} = nx Using Progression pr;
  Transfinite Line{l3,l5,l9,l10,l11,l12,l18,l21,l23,l25} = ny Using Progression pr;

  Transfinite Line{f1} = ny Using Progression pr;
  Transfinite Line{f2} = ny Using Progression pr;
  Transfinite Line{f3} = ny Using Progression pr;    
  Transfinite Line{f4} = ny Using Progression pr;
  Transfinite Line{f5} = ny Using Progression pr;
  Transfinite Line{f6} = ny Using Progression pr;
  Transfinite Line{f7} = nx Using Progression pr;
  Transfinite Line{f8} = ny Using Progression pr;  
  Transfinite Line{f9} = ny Using Progression pr;    
  Transfinite Line{f10} = ny Using Progression pr;
  Transfinite Line{f11} = ny Using Progression pr;
  Transfinite Line{f12} = ny Using Progression pr;
  Transfinite Line{f13} = ny Using Progression pr;
  Transfinite Line{f14} = ny Using Progression pr;


// Definição da superfície 

  Line Loop(1) = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25};
  Plane Surface(1) = {1};

  Line{f1} In Surface{1};
  Line{f2} In Surface{1};
  Line{f3} In Surface{1};
  Line{f4} In Surface{1};
  Line{f5} In Surface{1};
  Line{f6} In Surface{1};
  Line{f7} In Surface{1};
  Line{f8} In Surface{1};
  Line{f9} In Surface{1};
  Line{f10} In Surface{1};
  Line{f11} In Surface{1};
  Line{f12} In Surface{1};
  Line{f13} In Surface{1};
  Line{f14} In Surface{1};

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12};
  Physical Line("top") = {l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24};
  Physical Line("right") = {l13};
  Physical Line("left") = {l25};
  
  Physical Line("f1") = {f1};

  Physical Line("f2") = {f2};

  Physical Line("f3") = {f3};

  Physical Line("f4") = {f4};

  Physical Line("f5") = {f5};

  Physical Line("f6") = {f6};

  Physical Line("f7") = {f7};

  Physical Line("f8") = {f8};

  Physical Line("f9") = {f9};

  Physical Line("f10") = {f10};

  Physical Line("f11") = {f11};

  Physical Line("f12") = {f12};

  Physical Line("f13") = {f13};

  Physical Line("f14") = {f14};      
  
  Coherence Mesh;


