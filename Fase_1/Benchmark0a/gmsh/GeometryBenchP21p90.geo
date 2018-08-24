
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

  p5 = newp; Point(p5) = {16.2749023439828, 6.23486328125, 0, 1e+22};
  p6 = newp; Point(p6) = {13.8852539060172, 10, 0, 1e+22};
  p7 = newp; Point(p7) = {144.848144531017, 0, 0, 1e+22};
  p8 = newp; Point(p8) = {138.623535156017, 10, 0, 1e+22};
  p9 = newp; Point(p9) = {96.4799804689828, 0, 0, 1e+22};
  p10 = newp; Point(p10) = {91.3334960939828, 10, 0, 1e+22};
  p11 = newp; Point(p11) = {116.433105468983, 0, 0, 1e+22};
  p12 = newp; Point(p12) = {111.11328125, 10, 0, 1e+22};
  p13 = newp; Point(p13) = {175.595703125, 0, 0, 1e+22};
  p14 = newp; Point(p14) = {173.182617187966, 4.248046875, 0, 1e+22};
  p15 = newp; Point(p15) = {34.7060546879657, 0, 0, 1e+22};
  p16 = newp; Point(p16) = {32.5639648439828, 2.9326171875, 0, 1e+22};
  p17 = newp; Point(p17) = {101.584960937966, 0, 0, 1e+22};
  p18 = newp; Point(p18) = {99.2885742189828, 10, 0, 1e+22};
  p19 = newp; Point(p19) = {55.1469726560172, 0, 0, 1e+22};
  p20 = newp; Point(p20) = {51.3994140629657, 6.646484375, 0, 1e+22};
  p21 = newp; Point(p21) = {1.27636718796566, 0, 0, 1e+22};
  p22 = newp; Point(p22) = {0, 3.6279296875, 0, 1e+22};
  p23 = newp; Point(p23) = {45.599609375, 6.6376953125, 0, 1e+22};
  p24 = newp; Point(p24) = {44.9741210939828, 10, 0, 1e+22};
  p25 = newp; Point(p25) = {200, 9.69091796875, 0, 1e+22};
  p26 = newp; Point(p26) = {199.880859375, 10, 0, 1e+22};
  p27 = newp; Point(p27) = {153.142089843983, 0, 0, 1e+22};
  p28 = newp; Point(p28) = {150.48828125, 10, 0, 1e+22};
  p29 = newp; Point(p29) = {82.1235351560172, 10, 0, 1e+22};
  p30 = newp; Point(p30) = {83.0419921879657, 0, 0, 1e+22};
  p31 = newp; Point(p31) = {125.870117187034, 0, 0, 1e+22};
  p32 = newp; Point(p32) = {123.298828125, 10, 0, 1e+22};
  p33 = newp; Point(p33) = {172.943359375, 0, 0, 1e+22};
  p34 = newp; Point(p34) = {168.052734375, 10, 0, 1e+22};
  p35 = newp; Point(p35) = {66.0795898439828, 0, 0, 1e+22};
  p36 = newp; Point(p36) = {64.4150390629657, 10, 0, 1e+22};
  p37 = newp; Point(p37) = {13.9760742189828, 0, 0, 1e+22};
  p38 = newp; Point(p38) = {10.955078125, 10, 0, 1e+22};
  p39 = newp; Point(p39) = {4.24365234398283, 0, 0, 1e+22};
  p40 = newp; Point(p40) = {0, 6.3251953125, 0, 1e+22};
  p41 = newp; Point(p41) = {161.393554687966, 0, 0, 1e+22};
  p42 = newp; Point(p42) = {160.513671875, 6.15087890625, 0, 1e+22};
  p43 = newp; Point(p43) = {85.826171875, 0, 0, 1e+22};
  p44 = newp; Point(p44) = {82.869140625, 9.419921875, 0, 1e+22};
  p45 = newp; Point(p45) = {172.118652343983, 0, 0, 1e+22};
  p46 = newp; Point(p46) = {170.742675781017, 10, 0, 1e+22};
  p47 = newp; Point(p47) = {141.385742187034, 0, 0, 1e+22};
  p48 = newp; Point(p48) = {137.514160156017, 7.1572265625, 0, 1e+22};
  p49 = newp; Point(p49) = {145.985839843983, 4.21630859375, 0, 1e+22};
  p50 = newp; Point(p50) = {145.889648437034, 0, 0, 1e+22};

  p100 = newp;  Point(p100) = {171.7957815,2.346485134, 0, 1e+22};

// Fronteiras - Domínio Omega  
 
  l1 = newl; Line(l1) = {p1, p21};
  l2 = newl; Line(l2) = {p21, p39};
  l3 = newl; Line(l3) = {p39, p37};
  l4 = newl; Line(l4) = {p37, p15};
  l5 = newl; Line(l5) = {p15, p19};
  l6 = newl; Line(l6) = {p19, p35};
  l7 = newl; Line(l7) = {p35, p30};
  l8 = newl; Line(l8) = {p30, p43};
  l9 = newl; Line(l9) = {p43, p9};
  l10 = newl; Line(l10) = {p9, p17};
  l11 = newl; Line(l11) = {p17, p11};
  l12 = newl; Line(l12) = {p11, p31};
  l13 = newl; Line(l13) = {p31, p47};
  l14 = newl; Line(l14) = {p47, p7};
  l15 = newl; Line(l15) = {p7, p50};
  l16 = newl; Line(l16) = {p50, p27};
  l17 = newl; Line(l17) = {p27, p41};
  l18 = newl; Line(l18) = {p41, p45};
  l19 = newl; Line(l19) = {p45, p33};
  l20 = newl; Line(l20) = {p33, p13};
  l21 = newl; Line(l21) = {p13, p2};
  l22 = newl; Line(l22) = {p2, p25};
  l23 = newl; Line(l23) = {p25, p3};
  l24 = newl; Line(l24) = {p3, p26};
  l25 = newl; Line(l25) = {p26, p46};
  l26 = newl; Line(l26) = {p46, p34};
  l27 = newl; Line(l27) = {p34, p28};
  l28 = newl; Line(l28) = {p28, p8};
  l29 = newl; Line(l29) = {p8, p32};
  l30 = newl; Line(l30) = {p32, p12};
  l31 = newl; Line(l31) = {p12, p18};
  l32 = newl; Line(l32) = {p18, p10};
  l33 = newl; Line(l33) = {p10, p29};
  l34 = newl; Line(l34) = {p29, p36};
  l35 = newl; Line(l35) = {p36, p24};
  l36 = newl; Line(l36) = {p24, p6};
  l37 = newl; Line(l37) = {p6, p38};
  l38 = newl; Line(l38) = {p38, p4};
  l39 = newl; Line(l39) = {p4, p40};
  l40 = newl; Line(l40) = {p40, p22};
  l41 = newl; Line(l41) = {p22, p1};

//Fratura

  fshift = 41;

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
  f13 = newl; Line(f13) = {p29,p30};
  f14 = newl; Line(f14) = {p31,p32};
  
  f15a = newl; Line(f15a) = {p33,p100};
  f15b = newl; Line(f15b) = {p100,p34};

  f16 = newl; Line(f16) = {p35,p36};
  f17 = newl; Line(f17) = {p37,p38};
  f18 = newl; Line(f18) = {p39,p40};
  f19 = newl; Line(f19) = {p41,p42};
  f20 = newl; Line(f20) = {p43,p44};
  
  f21a = newl; Line(f21a) = {p45,p100};
  f21b = newl; Line(f21b) = {p100,p46};
  
  f22 = newl; Line(f22) = {p47,p48};
  f23 = newl; Line(f23) = {p49,p50};  

  Transfinite Line{l4,l5,l13,l21,l23,l25,l35,l36} = nx Using Progression pr;
  Transfinite Line{l1,l2,l3,l6,l7,l8,l9,l10,l11,l12,l14,l15,l16,l17,l18,l19,l20,l22,l24,l26,l27,l28,l29,l30,l31,l32,l33,l34,l37,l38,l39,l40,l41} = ny Using Progression pr;

  Transfinite Line{f1} = nf Using Progression pr;
  Transfinite Line{f2} = nf Using Progression pr;
  Transfinite Line{f3} = nf Using Progression pr;    
  Transfinite Line{f4} = nf Using Progression pr;
  Transfinite Line{f5} = nf Using Progression pr;
  Transfinite Line{f6} = nf Using Progression pr;
  Transfinite Line{f7} = nf Using Progression pr;
  Transfinite Line{f8} = nf Using Progression pr;  
  Transfinite Line{f9} = nf Using Progression pr;    
  Transfinite Line{f10} = nf Using Progression pr;
  Transfinite Line{f11} = nf Using Progression pr;
  Transfinite Line{f12} = nf Using Progression pr;
  Transfinite Line{f13} = nf Using Progression pr;
  Transfinite Line{f14} = nf Using Progression pr;
  Transfinite Line{f15a} = nf Using Progression pr;
  Transfinite Line{f15b} = nf Using Progression pr;
  Transfinite Line{f16} = nf Using Progression pr;
  Transfinite Line{f17} = nf Using Progression pr;
  Transfinite Line{f18} = nf Using Progression pr;
  Transfinite Line{f19} = nf Using Progression pr;
  Transfinite Line{f20} = nf Using Progression pr;
  Transfinite Line{f21a} = nf Using Progression pr;
  Transfinite Line{f21b} = nf Using Progression pr;
  Transfinite Line{f22} = nf Using Progression pr;  
  Transfinite Line{f23} = nf Using Progression pr;    


// Definição da superfície 

  Line Loop(1) = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38,l39,l40,l41};
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
  Line{f15a} In Surface{1};
  Line{f15b} In Surface{1};
  Line{f16} In Surface{1};
  Line{f17} In Surface{1};
  Line{f18} In Surface{1};
  Line{f19} In Surface{1};
  Line{f20} In Surface{1};
  Line{f21a} In Surface{1};
  Line{f21b} In Surface{1};
  Line{f22} In Surface{1};
  Line{f23} In Surface{1};

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21};
  Physical Line("right") = {l22,l23};
  Physical Line("top") = {l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38};
  Physical Line("left") = {l39,l40,l41};
  
  Physical Line("f1") = {f1};
//  Physical Point("P5") = {p5};
//  Physical Point("P6") = {p6};

  Physical Line("f2") = {f2};
//  Physical Point("P7") = {p7};
//  Physical Point("P8") = {p8};

  Physical Line("f3") = {f3};
//  Physical Point("P9") = {p9};
//  Physical Point("P10") = {p10};

  Physical Line("f4") = {f4};
//  Physical Point("P11") = {p11};
//  Physical Point("P12") = {p12};

  Physical Line("f5") = {f5};
//  Physical Point("P13") = {p13};
//  Physical Point("P14") = {p14};

  Physical Line("f6") = {f6};
//  Physical Point("P15") = {p15};
//  Physical Point("P16") = {p16};

  Physical Line("f7") = {f7};
//  Physical Point("P17") = {p17};
//  Physical Point("P18") = {p18};

  Physical Line("f8") = {f8};
//  Physical Point("P19") = {p19};
//  Physical Point("P20") = {p20};

  Physical Line("f9") = {f9};
//  Physical Point("P21") = {p21};
//  Physical Point("P22") = {p22};

  Physical Line("f10") = {f10};
//  Physical Point("P23") = {p23};
//  Physical Point("P24") = {p24};

  Physical Line("f11") = {f11};
//  Physical Point("P25") = {p25};
//  Physical Point("P26") = {p26};

  Physical Line("f12") = {f12};
//  Physical Point("P27") = {p27};
//  Physical Point("P28") = {p28};

  Physical Line("f13") = {f13};
//  Physical Point("P29") = {p29};
//  Physical Point("P30") = {p30};

  Physical Line("f14") = {f14};
//  Physical Point("P31") = {p31};
//  Physical Point("P32") = {p32};

  Physical Line("f15") = {f15a,f15b};      
//  Physical Point("P33") = {p33};
//  Physical Point("P100") = {p100};

//  Physical Line("f15") = {f15b};      
//  Physical Point("P100a") = {p100};
//  Physical Point("P34") = {p34};

  Physical Line("f16") = {f16};
//  Physical Point("P35") = {p35};
//  Physical Point("P36") = {p36};

  Physical Line("f17") = {f17};
//  Physical Point("P37") = {p37};
//  Physical Point("P38") = {p38};

  Physical Line("f18") = {f18};
//  Physical Point("P39") = {p39};
//  Physical Point("P40") = {p40};

  Physical Line("f19") = {f19};      
//  Physical Point("P41") = {p41};
//  Physical Point("P42") = {p42};

  Physical Line("f20") = {f20};
//  Physical Point("P43") = {p43};
//  Physical Point("P44") = {p44};

  Physical Line("f21") = {f21a,f21b};
//  Physical Point("P45") = {p45};
//  Physical Point("P100b") = {p100};

//  Physical Line("f21") = {f21b};
//  Physical Point("P100c") = {p100};
//  Physical Point("P46") = {p46};

  Physical Line("f22") = {f22};      
//  Physical Point("P47") = {p47};
//  Physical Point("P48") = {p48};

  Physical Line("f23") = {f23};      
//  Physical Point("P49") = {p49};
//  Physical Point("P50") = {p50};

  Coherence Mesh;



