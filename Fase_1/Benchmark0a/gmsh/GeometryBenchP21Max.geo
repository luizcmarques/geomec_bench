
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
nf = 8;
pr = 1;

// Coordenadas dos pontos

  //Domínio Omega
  p1 = newp; Point(p1) = {0, 0, 0, 1e+22};
  p2 = newp; Point(p2) = {L, 0, 0, 1e+22};
  p3 = newp; Point(p3) = {L, h, 0, 1e+22};
  p4 = newp; Point(p4) = {0, h, 0, 1e+22};

  //Fraturas
  
  dshift=4;

  p5 = newp; Point(p5) = {26.70068359375, 10, 0, 1e+22};
  p6 = newp; Point(p6) = {29.08984375, 6.23486328125, 0, 1e+22};
  p7 = newp; Point(p7) = {93.11572265625, 10, 0, 1e+22}; 
  p8 = newp; Point(p8) = {94.50732421875, 0, 0, 1e+22}; 
  p9 = newp; Point(p9) = {122.5947265625, 10, 0, 1e+22}; 
  p10 = newp; Point(p10) = {127.7412109375, 0, 0, 1e+22}; 
  p11 = newp; Point(p11) = {193.3662109375, 9.9541015625, 0, 1e+22}; 
  p12 = newp; Point(p12) = {200, 1.23388671875, 0, 1e+22}; 
  p13 = newp; Point(p13) = {154.40966796875, 10, 0, 1e+22}; 
  p14 = newp; Point(p14) = {160.11865234375, 0, 0, 1e+22}; 
  p15 = newp; Point(p15) = {56.56591796875, 10, 0, 1e+22}; 
  p16 = newp; Point(p16) = {60.287109375, 0, 0, 1e+22}; 
  p17 = newp; Point(p17) = {120.54345703125, 10, 0, 1e+22}; 
  p18 = newp; Point(p18) = {125.86328125,  0, 0, 1e+22}; 
  p19 = newp; Point(p19) = {190.78125, 4.248046875, 0, 1e+22}; 
  p20 = newp; Point(p20) = {193.1943359375, 0, 0, 1e+22}; 
  p21 = newp; Point(p21) = {51.29345703125, 2.9326171875, 0, 1e+22}; 
  p22 = newp; Point(p22) = {53.435546875, 0, 0, 1e+22}; 

  p23 = newp; Point(p23) = {120.6572265625, 10, 0, 1e+22}; 
  p24 = newp; Point(p24) = {122.95361328125, 0, 0, 1e+22}; 
  p25 = newp; Point(p25) = {138.86767578125, 6.60302734375, 0, 1e+22}; 
  p26 = newp; Point(p26) = {141.1396484375, 0, 0, 1e+22}; 
  p27 = newp; Point(p27) = {68.560546875, 6.646484375, 0, 1e+22}; 
  p28 = newp; Point(p28) = {72.30810546875, 0, 0, 1e+22}; 
  p29 = newp; Point(p29) = {7.595703125, 10, 0, 1e+22}; 
  p30 = newp; Point(p30) = {11.11376953125, 0, 0, 1e+22}; 
  p31 = newp; Point(p31) = {189.45458984375, 10, 0, 1e+22}; 
  p32 = newp; Point(p32) = {193.30615234375, 0, 0, 1e+22}; 
  p33 = newp; Point(p33) = {197.37255859375, 10, 0, 1e+22}; 
  p34 = newp; Point(p34) = {199.87548828125, 0, 0, 1e+22}; 

  p35 = newp; Point(p35) = {159.3076171875, 10, 0, 1e+22}; 
  p36 = newp; Point(p36) = {159.8720703125, 6.67578125, 0, 1e+22}; 
  p37 = newp; Point(p37) = {66.70068359375, 10, 0, 1e+22}; 
  p38 = newp; Point(p38) = {67.61865234375,  0, 0, 1e+22}; 
  p39 = newp; Point(p39) = {117.12451171875, 10, 0, 1e+22}; 
  p40 = newp; Point(p40) = {119.69580078125, 0, 0, 1e+22}; 
  p41 = newp; Point(p41) = {150.11669921875, 10, 0, 1e+22}; 
  p42 = newp; Point(p42) = {155.00732421875, 0, 0, 1e+22}; 
  p43 = newp; Point(p43) = {130.556640625, 10, 0, 1e+22}; 
  p44 = newp; Point(p44) = {133.2509765625,  0, 0, 1e+22}; 
  p45 = newp; Point(p45) = {70.60595703125, 10, 0, 1e+22}; 
  p46 = newp; Point(p46) = {72.8486328125, 0, 0, 1e+22}; 
  p47 = newp; Point(p47) = {52.7529296875, 10, 0, 1e+22}; 
  p48 = newp; Point(p48) = {54.4169921875, 0, 0, 1e+22}; 
  p49 = newp; Point(p49) = {0.5771484375, 10, 0, 1e+22}; 
  p50 = newp; Point(p50) = {4.5419921875, 0, 0, 1e+22}; 
  
  p51 = newp; Point(p51) = {4.32177734375, 10, 0, 1e+22}; 
  p52 = newp; Point(p52) = {9.75439453125, 0, 0, 1e+22}; 
  p53 = newp; Point(p53) = {5.6865234375, 10, 0, 1e+22}; 
  p54 = newp; Point(p54) = {8.70703125, 0, 0, 1e+22}; 

  p55 = newp; Point(p55) = {63.36669921875, 10, 0, 1e+22}; 
  p56 = newp; Point(p56) = {68.806640625, 0, 0, 1e+22}; 
  p57 = newp; Point(p57) = {147.62890625, 6.15087890625, 0, 1e+22}; 
  p58 = newp; Point(p58) = {148.50830078125, 0, 0, 1e+22}; 
  p59 = newp; Point(p59) = {76.34716796875, 9.419921875, 0, 1e+22}; 
  p60 = newp; Point(p60) = {79.30419921875, 0, 0, 1e+22}; 
  p61 = newp; Point(p61) = {119.36376953125, 7.1572265625, 0, 1e+22}; 
  p62 = newp; Point(p62) = {123.23486328125, 0, 0, 1e+22};
  p63 = newp; Point(p63) = {139.31201171875, 0, 0, 1e+22}; 
  p64 = newp; Point(p64) = {139.408203125, 4.21630859375, 0, 1e+22}; 

  p100 = newp; Point(p100) = {7.395494962993421, 4.342105263157896, 0, 1e+22}; //Intersection p65
  p101 = newp; Point(p101) = {67.37748854191973, 2.627146097, 0, 1e+22}; //Intersection
  p102 = newp; Point(p102) = {122.74609271220551,0.9036830223, 0, 1e+22}; //Intersection
  p103 = newp; Point(p103) = {120.7436378,9.623708010335909, 0, 1e+22}; //Intersection
  p104 = newp; Point(p104) = {199.35416716586613,2.08284362915751, 0, 1e+22}; //Intersection

// Fronteiras - Domínio Omega  

  l1 = newl; Line(l1) = {p1,p50}; 
  l2 = newl; Line(l2) = {p50,p54};
  l3 = newl; Line(l3) = {p54,p52};
  l4 = newl; Line(l4) = {p52,p30};
  l5 = newl; Line(l5) = {p30,p22};
  l6 = newl; Line(l6) = {p22,p48};
  l7 = newl; Line(l7) = {p48,p16};
  l8 = newl; Line(l8) = {p16,p38};
  l9 = newl; Line(l9) = {p38,p56};
  l10 = newl; Line(l10) = {p56,p28};
  l11 = newl; Line(l11) = {p28,p46};
  l12 = newl; Line(l12) = {p46,p60};
  l13 = newl; Line(l13) = {p60,p8};
  l14 = newl; Line(l14) = {p8,p40};
  l15 = newl; Line(l15) = {p40,p24};
  l16 = newl; Line(l16) = {p24,p62};
  l17 = newl; Line(l17) = {p62,p18};
  l18 = newl; Line(l18) = {p18,p10};
  l19 = newl; Line(l19) = {p10,p44};
  l20 = newl; Line(l20) = {p44,p63};
  l21 = newl; Line(l21) = {p63,p26};
  l22 = newl; Line(l22) = {p26,p58};
  l23 = newl; Line(l23) = {p58,p42};
  l24 = newl; Line(l24) = {p42,p14};
  l25 = newl; Line(l25) = {p14,p20};
  l26 = newl; Line(l26) = {p32,p34};
  l27 = newl; Line(l27) = {p20,p32};
  l28 = newl; Line(l28) = {p34,p2};

  l29 = newl; Line(l29) = {p2,p12};
  l30 = newl; Line(l30) = {p12,p3};
 
  l31 = newl; Line(l31) = {p3,p33};
  l32 = newl; Line(l32) = {p33,p31};
  l33 = newl; Line(l33) = {p31,p35};
  l34 = newl; Line(l34) = {p35,p13};
  l35 = newl; Line(l35) = {p13,p41};
  l36 = newl; Line(l36) = {p41,p43};
  l37 = newl; Line(l37) = {p43,p9};
  l38 = newl; Line(l38) = {p9,p23}; 
  l39 = newl; Line(l39) = {p23,p17};
  l40 = newl; Line(l40) = {p17,p39};

  l41 = newl; Line(l41) = {p39, p7}; 
  l42 = newl; Line(l42) = {p7, p45};
  l43 = newl; Line(l43) = {p45, p37};
  l44 = newl; Line(l44) = {p37, p55};
  l45 = newl; Line(l45) = {p55, p15};
  l46 = newl; Line(l46) = {p15, p47};
  l47 = newl; Line(l47) = {p47, p5};
  l48 = newl; Line(l48) = {p5, p29};
  l49 = newl; Line(l49) = {p29, p53};
  l50 = newl; Line(l50) = {p53, p51};
  l51 = newl; Line(l51) = {p51, p49};
  l52 = newl; Line(l52) = {p49, p4};
  
  l53 = newl; Line(l53) = {p4,p1};

//Fratura

  fshift = 53;

  f1 = newl; Line(f1) = {p5,p6};
  f2 = newl; Line(f2) = {p7,p8};
  f3 = newl; Line(f3) = {p9,p10};
  
  f4a = newl; Line(f4a) = {p11,p104};
  f4b = newl; Line(f4b) = {p104,p12};

  f5 = newl; Line(f5) = {p13,p14};
  f6 = newl; Line(f6) = {p15,p16};
  
  f7a = newl; Line(f7a) = {p17,p103};
  f7b = newl; Line(f7b) = {p103,p18};
  
  f8 = newl; Line(f8) = {p19,p20};
  f9 = newl; Line(f9) = {p21,p22};

  f10a = newl; Line(f10a) = {p23,p103};
  f10b = newl; Line(f10b) = {p102,p103};
  f10c = newl; Line(f10c) = {p102,p24};

  f11 = newl; Line(f11) = {p25,p26};
  f12 = newl; Line(f12) = {p27,p28};
  f13 = newl; Line(f13) = {p29,p30};
  f14 = newl; Line(f14) = {p31,p32};
  
  f15a = newl; Line(f15a) = {p33,p104};
  f15b = newl; Line(f15b) = {p104,p34};

  f16 = newl; Line(f16) = {p35,p36};
  
  f17a = newl; Line(f17a) = {p37,p101};
  f17b = newl; Line(f17b) = {p101,p38};

  f18 = newl; Line(f18) = {p39,p40};
  f19 = newl; Line(f19) = {p41,p42};
  f20 = newl; Line(f20) = {p43,p44};
  f21 = newl; Line(f21) = {p45,p46};
  f22 = newl; Line(f22) = {p47,p48};
  f23 = newl; Line(f23) = {p49,p50};
  
  f24a = newl; Line(f24a) = {p51,p100};
  f24b = newl; Line(f24b) = {p100,p52};

  f25a = newl; Line(f25a) = {p53,p100};
  f25b = newl; Line(f25b) = {p100,p54};
  
  f26a = newl; Line(f26a) = {p55,p101};
  f26b = newl; Line(f26b) = {p101,p56};

  f27 = newl; Line(f27) = {p57,p58};
  f28 = newl; Line(f28) = {p59,p60};
  
  f29a = newl; Line(f29a) = {p61,p102};
  f29b = newl; Line(f29b) = {p102,p62};

  f30 = newl; Line(f30) = {p63,p64};


  Transfinite Line{l5,l13,l14,l22,l24,l25,l33,l36,l41,l42,l47,l48,l53} = nx Using Progression pr;
  Transfinite Line{l1,l2,l3,l4,l6,l7,l8,l9,l10,l11,l12,l15,l16,l17,l18,l19,l20,l21,l23,l26,l29,l30,l31,l32,l34,l35,l37,l38,l39,l40,l43,l44,l45,l46,l49,l50,l51,l52} = ny Using Progression pr;

  Transfinite Line{f1} = nf Using Progression pr;
  Transfinite Line{f2} = nf Using Progression pr;
  Transfinite Line{f3} = nf Using Progression pr;    
  Transfinite Line{f4a} = nf Using Progression pr;
  Transfinite Line{f4b} = nf Using Progression pr;
  Transfinite Line{f5} = nf Using Progression pr;
  Transfinite Line{f6} = nf Using Progression pr;
  Transfinite Line{f7a} = nf Using Progression pr;
  Transfinite Line{f7b} = nf Using Progression pr;
  Transfinite Line{f8} = nf Using Progression pr;  
  Transfinite Line{f9} = nf Using Progression pr;    
  Transfinite Line{f10a} = nf Using Progression pr;
  Transfinite Line{f10b} = nf Using Progression pr;
  Transfinite Line{f10c} = nf Using Progression pr;
  Transfinite Line{f11} = nf Using Progression pr;
  Transfinite Line{f12} = nf Using Progression pr;
  Transfinite Line{f13} = nf Using Progression pr;
  Transfinite Line{f14} = nf Using Progression pr;
  Transfinite Line{f15a} = nf Using Progression pr;
  Transfinite Line{f15b} = nf Using Progression pr;
  Transfinite Line{f16} = nf Using Progression pr;
  Transfinite Line{f17a} = nf Using Progression pr;
  Transfinite Line{f17b} = nf Using Progression pr;    
  Transfinite Line{f18} = nf Using Progression pr;
  Transfinite Line{f19} = nf Using Progression pr;
  Transfinite Line{f20} = nf Using Progression pr;
  Transfinite Line{f21} = nf Using Progression pr;
  Transfinite Line{f22} = nf Using Progression pr;  
  Transfinite Line{f23} = nf Using Progression pr;    
  Transfinite Line{f24a} = nf Using Progression pr;
  Transfinite Line{f24b} = nf Using Progression pr;
  Transfinite Line{f25a} = nf Using Progression pr;
  Transfinite Line{f25b} = nf Using Progression pr;
  Transfinite Line{f26a} = nf Using Progression pr;
  Transfinite Line{f26b} = nf Using Progression pr;
  Transfinite Line{f27} = nf Using Progression pr;
  Transfinite Line{f28} = nf Using Progression pr;
  Transfinite Line{f29a} = nf Using Progression pr;
  Transfinite Line{f29b} = nf Using Progression pr;
  Transfinite Line{f30} = nf Using Progression pr;


// Definição da superfície 

  Line Loop(1) = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38,l39,l40,l41,l42,l43,l44,l45,l46,l47,l48,l49,l50,l51,l52,l53};
  Plane Surface(1) = {1};

  Line{f1} In Surface{1};
  Line{f2} In Surface{1};
  Line{f3} In Surface{1};
  Line{f4a} In Surface{1};
  Line{f4b} In Surface{1};
  Line{f5} In Surface{1};
  Line{f6} In Surface{1};
  Line{f7a} In Surface{1};
  Line{f7b} In Surface{1};
  Line{f8} In Surface{1};
  Line{f9} In Surface{1};
  Line{f10a} In Surface{1};
  Line{f10b} In Surface{1};
  Line{f10c} In Surface{1};
  Line{f11} In Surface{1};
  Line{f12} In Surface{1};
  Line{f13} In Surface{1};
  Line{f14} In Surface{1};
  Line{f15a} In Surface{1};
  Line{f15b} In Surface{1};
  Line{f16} In Surface{1};
  Line{f17a} In Surface{1};
  Line{f17b} In Surface{1};
  Line{f18} In Surface{1};
  Line{f19} In Surface{1};
  Line{f20} In Surface{1};
  Line{f21} In Surface{1};
  Line{f22} In Surface{1};
  Line{f23} In Surface{1};
  Line{f24a} In Surface{1};
  Line{f24b} In Surface{1};
  Line{f25a} In Surface{1};
  Line{f25b} In Surface{1};
  Line{f26a} In Surface{1};
  Line{f26b} In Surface{1};
  Line{f27} In Surface{1};
  Line{f28} In Surface{1};
  Line{f29a} In Surface{1};
  Line{f29b} In Surface{1};
  Line{f30} In Surface{1};

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26,l27,l28};
  Physical Line("top") = {l31,l32,l33,l34,l35,l36,l37,l38,l39,l40,l41,l42,l43,l44,l45,l46,l47,l48,l49,l50,l51,l52};
  Physical Line("right") = {l29,l30};
  Physical Line("left") = {l53};
  
  Physical Line("f1") = {f1};

  Physical Line("f2") = {f2};

  Physical Line("f3") = {f3};

  Physical Line("f4") = {f4a,f4b};

  Physical Line("f5") = {f5};

  Physical Line("f6") = {f6};

  Physical Line("f7") = {f7a,f7b};

  Physical Line("f8") = {f8};

  Physical Line("f9") = {f9};

  Physical Line("f10a") = {f10a,f10b,f10c};

  Physical Line("f11") = {f11};

  Physical Line("f12") = {f12};

  Physical Line("f13") = {f13};

  Physical Line("f14") = {f14};

  Physical Line("f15") = {f15a,f15b};      

  Physical Line("f16") = {f16};

  Physical Line("f17") = {f17a,f17b};

  Physical Line("f18") = {f18};

  Physical Line("f19") = {f19};      

  Physical Line("f20") = {f20};

  Physical Line("f21") = {f21};

  Physical Line("f22") = {f22};      

  Physical Line("f23") = {f23};      

  Physical Line("f24") = {f24a,f24b};      

  Physical Line("f25") = {f25a,f25b};      

  Physical Line("f26") = {f26a,f26b};      

  Physical Line("f27") = {f27};      

  Physical Line("f28") = {f28};      

  Physical Line("f29") = {f29a,f29b};      

  Physical Line("f30") = {f30};      

  Coherence Mesh;


