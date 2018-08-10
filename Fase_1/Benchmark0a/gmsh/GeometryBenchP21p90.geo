
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
  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {L, 0, 0, 1e+22};
  Point(3) = {L, h, 0, 1e+22};
  Point(4) = {0, h, 0, 1e+22};

  //Fraturas
  
  dshift=4;

  Point(1+dshift) = {16.2749023439828, 6.23486328125, 0, 1e+22};
  Point(2+dshift) = {13.8852539060172, 10, 0, 1e+22};
  Point(3+dshift) = {144.848144531017, 0, 0, 1e+22};
  Point(4+dshift) = {138.623535156017, 10, 0, 1e+22};
  Point(5+dshift) = {96.4799804689828, 0, 0, 1e+22};
  Point(6+dshift) = {91.3334960939828, 10, 0, 1e+22};
  Point(7+dshift) = {116.433105468983, 0, 0, 1e+22};
  Point(8+dshift) = {111.11328125, 10, 0, 1e+22};
  Point(9+dshift) = {175.595703125, 0, 0, 1e+22};
  Point(10+dshift) = {173.182617187966, 4.248046875, 0, 1e+22};
  Point(11+dshift) = {34.7060546879657, 0, 0, 1e+22};
  Point(12+dshift) = {32.5639648439828, 2.9326171875, 0, 1e+22};
  Point(13+dshift) = {101.584960937966, 0, 0, 1e+22};
  Point(14+dshift) = {99.2885742189828, 10, 0, 1e+22};
  Point(15+dshift) = {55.1469726560172, 0, 0, 1e+22};
  Point(16+dshift) = {51.3994140629657, 6.646484375, 0, 1e+22};
  Point(17+dshift) = {1.27636718796566, 0, 0, 1e+22};
  Point(18+dshift) = {0, 3.6279296875, 0, 1e+22};
  Point(19+dshift) = {45.599609375, 6.6376953125, 0, 1e+22};
  Point(20+dshift) = {44.9741210939828, 10, 0, 1e+22};
  Point(21+dshift) = {200, 9.69091796875, 0, 1e+22};
  Point(22+dshift) = {199.880859375, 10, 0, 1e+22};
  Point(23+dshift) = {153.142089843983, 0, 0, 1e+22};
  Point(24+dshift) = {150.48828125, 10, 0, 1e+22};
  Point(25+dshift) = {82.1235351560172, 10, 0, 1e+22};
  Point(26+dshift) = {83.0419921879657, 0, 0, 1e+22};
  Point(27+dshift) = {125.870117187034, 0, 0, 1e+22};
  Point(28+dshift) = {123.298828125, 10, 0, 1e+22};
  Point(29+dshift) = {172.943359375, 0, 0, 1e+22};
  Point(30+dshift) = {168.052734375, 10, 0, 1e+22};
  Point(31+dshift) = {66.0795898439828, 0, 0, 1e+22};
  Point(32+dshift) = {64.4150390629657, 10, 0, 1e+22};
  Point(33+dshift) = {13.9760742189828, 0, 0, 1e+22};
  Point(34+dshift) = {10.955078125, 10, 0, 1e+22};
  Point(35+dshift) = {4.24365234398283, 0, 0, 1e+22};
  Point(36+dshift) = {0, 6.3251953125, 0, 1e+22};
  Point(37+dshift) = {161.393554687966, 0, 0, 1e+22};
  Point(38+dshift) = {160.513671875, 6.15087890625, 0, 1e+22};
  Point(39+dshift) = {85.826171875, 0, 0, 1e+22};
  Point(40+dshift) = {82.869140625, 9.419921875, 0, 1e+22};
  Point(41+dshift) = {172.118652343983, 0, 0, 1e+22};
  Point(42+dshift) = {170.742675781017, 10, 0, 1e+22};
  Point(43+dshift) = {141.385742187034, 0, 0, 1e+22};
  Point(44+dshift) = {137.514160156017, 7.1572265625, 0, 1e+22};
  Point(45+dshift) = {145.985839843983, 4.21630859375, 0, 1e+22};
  Point(46+dshift) = {145.889648437034, 0, 0, 1e+22};

  Point(100) = {171.7957815,2.346485134, 0, 1e+22};

// Fronteiras - Domínio Omega  
 
  Line(1) = {1, 21};
  Line(2) = {21, 39};
  Line(3) = {39, 37};
  Line(4) = {37, 15};
  Line(5) = {15, 19};
  Line(6) = {19, 35};
  Line(7) = {35, 30};
  Line(8) = {30, 43};
  Line(9) = {43, 9};
  Line(10) = {9, 17};
  Line(11) = {17, 11};
  Line(12) = {11, 31};
  Line(13) = {31, 47};
  Line(14) = {47, 7};
  Line(15) = {7, 50};
  Line(16) = {50, 27};
  Line(17) = {27, 41};
  Line(18) = {41, 45};
  Line(19) = {45, 33};
  Line(20) = {33, 13};
  Line(21) = {13, 2};
  Line(22) = {2, 25};
  Line(23) = {25, 3};
  Line(24) = {3, 26};
  Line(25) = {26, 46};
  Line(26) = {46, 34};
  Line(27) = {34, 28};
  Line(28) = {28, 8};
  Line(29) = {8, 32};
  Line(30) = {32, 12};
  Line(31) = {12, 18};
  Line(32) = {18, 10};
  Line(33) = {10, 29};
  Line(34) = {29, 36};
  Line(35) = {36, 24};
  Line(36) = {24, 6};
  Line(37) = {6, 38};
  Line(38) = {38, 4};
  Line(39) = {4, 40};
  Line(40) = {40, 22};
  Line(41) = {22, 1};

//Fratura

  fshift = 41;

  Line(1+fshift) = {1+dshift,2+dshift};
  Line(2+fshift) = {3+dshift,4+dshift};
  Line(3+fshift) = {5+dshift,6+dshift};
  Line(4+fshift) = {7+dshift,8+dshift};
  Line(5+fshift) = {9+dshift,10+dshift};
  Line(6+fshift) = {11+dshift,12+dshift};
  Line(7+fshift) = {13+dshift,14+dshift};
  Line(8+fshift) = {15+dshift,16+dshift};
  Line(9+fshift) = {17+dshift,18+dshift};
  Line(10+fshift) = {19+dshift,20+dshift};
  Line(11+fshift) = {21+dshift,22+dshift};
  Line(12+fshift) = {23+dshift,24+dshift};
  Line(13+fshift) = {25+dshift,26+dshift};
  Line(14+fshift) = {27+dshift,28+dshift};
  
  Line(15+fshift+100) = {29+dshift,100};
  Line(15+fshift+101) = {100,30+dshift};

  Line(16+fshift) = {31+dshift,32+dshift};
  Line(17+fshift) = {33+dshift,34+dshift};
  Line(18+fshift) = {35+dshift,36+dshift};
  Line(19+fshift) = {37+dshift,38+dshift};
  Line(20+fshift) = {39+dshift,40+dshift};
  
  Line(21+fshift+102) = {41+dshift,100};
  Line(21+fshift+103) = {100,42+dshift};
  
  Line(22+fshift) = {43+dshift,44+dshift};
  Line(23+fshift) = {45+dshift,46+dshift};  

  Transfinite Line{4,5,13,21,23,25,35,36} = nx Using Progression pr;
  Transfinite Line{1,2,3,6,7,8,9,10,11,12,14,15,16,17,18,19,20,22,24,26,27,28,29,30,31,32,33,34,37,38,39,40,41} = ny Using Progression pr;

  Transfinite Line{1+fshift} = ny Using Progression pr;
  Transfinite Line{2+fshift} = ny Using Progression pr;
  Transfinite Line{3+fshift} = ny Using Progression pr;    
  Transfinite Line{4+fshift} = ny Using Progression pr;
  Transfinite Line{5+fshift} = ny Using Progression pr;
  Transfinite Line{6+fshift} = ny Using Progression pr;
  Transfinite Line{7+fshift} = ny Using Progression pr;
  Transfinite Line{8+fshift} = ny Using Progression pr;  
  Transfinite Line{9+fshift} = ny Using Progression pr;    
  Transfinite Line{10+fshift} = ny Using Progression pr;
  Transfinite Line{11+fshift} = ny Using Progression pr;
  Transfinite Line{12+fshift} = ny Using Progression pr;
  Transfinite Line{13+fshift} = ny Using Progression pr;
  Transfinite Line{14+fshift} = ny Using Progression pr;
  Transfinite Line{15+fshift+100} = ny Using Progression pr;
  Transfinite Line{15+fshift+101} = ny Using Progression pr;
  Transfinite Line{16+fshift} = ny Using Progression pr;
  Transfinite Line{17+fshift} = ny Using Progression pr;
  Transfinite Line{18+fshift} = ny Using Progression pr;
  Transfinite Line{19+fshift} = ny Using Progression pr;
  Transfinite Line{20+fshift} = ny Using Progression pr;
  Transfinite Line{21+fshift+102} = ny Using Progression pr;
  Transfinite Line{21+fshift+103} = ny Using Progression pr;
  Transfinite Line{22+fshift} = ny Using Progression pr;  
  Transfinite Line{23+fshift} = ny Using Progression pr;    


// Definição da superfície 

  Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41};
  Plane Surface(1) = {1};

  Line{1+fshift} In Surface{1};
  Line{2+fshift} In Surface{1};
  Line{3+fshift} In Surface{1};
  Line{4+fshift} In Surface{1};
  Line{5+fshift} In Surface{1};
  Line{6+fshift} In Surface{1};
  Line{7+fshift} In Surface{1};
  Line{8+fshift} In Surface{1};
  Line{9+fshift} In Surface{1};
  Line{10+fshift} In Surface{1};
  Line{11+fshift} In Surface{1};
  Line{12+fshift} In Surface{1};
  Line{13+fshift} In Surface{1};
  Line{14+fshift} In Surface{1};
  Line{15+fshift+100} In Surface{1};
  Line{15+fshift+101} In Surface{1};
  Line{16+fshift} In Surface{1};
  Line{17+fshift} In Surface{1};
  Line{18+fshift} In Surface{1};
  Line{19+fshift} In Surface{1};
  Line{20+fshift} In Surface{1};
  Line{21+fshift+102} In Surface{1};
  Line{21+fshift+103} In Surface{1};
  Line{22+fshift} In Surface{1};
  Line{23+fshift} In Surface{1};

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
  Physical Line("top") = {25,26,27,28,29,30,31,32,33,34,35,36,37,38};
  Physical Line("right") = {22,23};
  Physical Line("left") = {39,40,41};
  
  Physical Line("f1") = {1+fshift};
  Physical Point("P5") = {5};
  Physical Point("P6") = {6};

  Physical Line("f2") = {2+fshift};
  Physical Point("P7") = {7};
  Physical Point("P8") = {8};

  Physical Line("f3") = {3+fshift};
  Physical Point("P9") = {9};
  Physical Point("P10") = {10};

  Physical Line("f4") = {4+fshift};
  Physical Point("P11") = {11};
  Physical Point("P12") = {12};

  Physical Line("f5") = {5+fshift};
  Physical Point("P13") = {13};
  Physical Point("P14") = {14};

  Physical Line("f6") = {6+fshift};
  Physical Point("P15") = {15};
  Physical Point("P16") = {16};

  Physical Line("f7") = {7+fshift};
  Physical Point("P17") = {17};
  Physical Point("P18") = {18};

  Physical Line("f8") = {8+fshift};
  Physical Point("P19") = {19};
  Physical Point("P20") = {20};

  Physical Line("f9") = {9+fshift};
  Physical Point("P21") = {21};
  Physical Point("P22") = {22};

  Physical Line("f10") = {10+fshift};
  Physical Point("P23") = {23};
  Physical Point("P24") = {24};

  Physical Line("f11") = {11+fshift};
  Physical Point("P25") = {25};
  Physical Point("P26") = {26};

  Physical Line("f12") = {12+fshift};
  Physical Point("P27") = {27};
  Physical Point("P28") = {28};

  Physical Line("f13") = {13+fshift};
  Physical Point("P29") = {29};
  Physical Point("P30") = {30};

  Physical Line("f14") = {14+fshift};
  Physical Point("P31") = {31};
  Physical Point("P32") = {32};

  Physical Line("f15") = {15+fshift+100};      
  Physical Point("P33") = {33};
  Physical Point("P100") = {100};

  Physical Line("f15b") = {15+fshift+101};      
  Physical Point("P100a") = {100};
  Physical Point("P34") = {34};

  Physical Line("f16") = {16+fshift};
  Physical Point("P35") = {35};
  Physical Point("P36") = {36};

  Physical Line("f17") = {17+fshift};
  Physical Point("P37") = {37};
  Physical Point("P38") = {38};

  Physical Line("f18") = {18+fshift};
  Physical Point("P39") = {39};
  Physical Point("P40") = {40};

  Physical Line("f19") = {19+fshift};      
  Physical Point("P41") = {41};
  Physical Point("P42") = {42};

  Physical Line("f20") = {20+fshift};
  Physical Point("P43") = {43};
  Physical Point("P44") = {44};

  Physical Line("f21") = {21+fshift+102};
  Physical Point("P45") = {45};
  Physical Point("P100b") = {100};

  Physical Line("f21b") = {21+fshift+103};
  Physical Point("P100c") = {100};
  Physical Point("P46") = {46};

  Physical Line("f22") = {22+fshift};      
  Physical Point("P47") = {47};
  Physical Point("P48") = {48};

  Physical Line("f23") = {23+fshift};      
  Physical Point("P49") = {49};
  Physical Point("P50") = {50};

  Coherence Mesh;



