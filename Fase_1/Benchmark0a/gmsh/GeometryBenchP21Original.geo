

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

  Point(5) = {59.56005859375, 0, 0, 1e+22}; 
  Point(6) = {57.26318359375, 10, 0, 1e+22};

  Point(7) = {76.92333984375, 5.39453125, 0, 1e+22};
  Point(8) = {76.5068359375, 10, 0, 1e+22};

  Point(9) = {89.01513671875, 0, 0, 1e+22};
  Point(10) = {87.35107421875, 10, 0, 1e+22};

  Point(11) = {97.42724609375,  0, 0, 1e+22};
  Point(12) = {90.37890625, 10, 0, 1e+22};

  Point(13) = {102.81201171875, 10, 0, 1e+22};
  Point(14) = {102.45849609375, 2.873046875, 0, 1e+22};

  Point(15) = {127.0439453125,  0, 0, 1e+22};
  Point(16) = {121.64697265625, 10, 0, 1e+22};

  Point(17) = {191.47509765625, 0, 0, 1e+22};
  Point(18) = {188.8212890625,  10, 0, 1e+22};

  Point(19) = {157.5068359375,  1.970703125, 0, 1e+22};
  Point(20) = {157.3876953125,  0, 0, 1e+22};

  Point(21) = {177.08203125,  0, 0, 1e+22};
  Point(22) = {173.21044921875, 7.1572265625, 0, 1e+22};

  Point(23) = {184.17529296875, 0, 0, 1e+22};
  Point(24) = {182.79931640625, 10, 0, 1e+22};

  Point(25) = {131.27392578125, 0, 0, 1e+22};
  Point(26) = {125.04931640625, 10, 0, 1e+22};

  Point(27) = {192.71142578125, 0, 0, 1e+22};
  Point(28) = {192.60400390625, 1.470703125, 0, 1e+22};

//  Point(29) = {87.95263671875,  6.6376953125, 0, 1e+22};
//  Point(30) = {87.3271484375, 10, 0, 1e+22};

  Point(31) = {193.7646484375,  3.07763671875, 0, 1e+22};
  Point(32) = {192.798828125, 10, 0, 1e+22};

  Point(33) = {195.78076171875, 0, 0, 1e+22};
  Point(34) = {195.13818359375, 1.36083984375, 0, 1e+22};

// Fronteiras

  //Domínio Omega  
  Line(1) = {1,5};
  Line(2) = {5,9};
  Line(3) = {9,11};
  Line(4) = {11,15};
  Line(5) = {15,25};  
  Line(6) = {25,20};
  Line(7) = {20,21};
  Line(8) = {21,23};
  Line(9) = {23,17};
  Line(10) = {17,27};
  Line(11) = {27,33};
  Line(12) = {33,2};
  Line(13) = {2,3};
  Line(14) = {3,32};
  Line(15) = {32,18};
  Line(16) = {18,24};
  Line(17) = {24,26};
  Line(18) = {26,16};
  Line(19) = {16,13};
  Line(20) = {13,12};
  Line(21) = {12,10};
  Line(22) = {10,8};
  Line(23) = {8,6};
  Line(24) = {6,4};
  Line(25) = {4,1};

  //Fratura

  fshift = 25;

  Line(1+fshift) = {5, 6};
  Line(2+fshift) = {7, 8};
  Line(3+fshift) = {9, 10};
  Line(4+fshift) = {11, 12};
  Line(5+fshift) = {13, 14};
  Line(6+fshift) = {15, 16};
  Line(7+fshift) = {17, 18};
  Line(8+fshift) = {19, 20};
  Line(9+fshift) = {21, 22};
  Line(10+fshift) = {23, 24};
  Line(11+fshift) = {25, 26};
  Line(12+fshift) = {27, 28};
//  Line() = {29, 30};
  Line(13+fshift) = {31, 32};
  Line(14+fshift) = {33, 34};

  Transfinite Line{1,2,4,6,7,8,13,14,15,16,17,19,20,22,24} = nx Using Progression pr;
  Transfinite Line{3,5,9,10,11,12,18,21,23,25} = ny Using Progression pr;

  Transfinite Line{1+fshift} = ny Using Progression pr;
  Transfinite Line{2+fshift} = ny Using Progression pr;
  Transfinite Line{3+fshift} = ny Using Progression pr;    
  Transfinite Line{4+fshift} = ny Using Progression pr;
  Transfinite Line{5+fshift} = ny Using Progression pr;
  Transfinite Line{6+fshift} = ny Using Progression pr;
  Transfinite Line{7+fshift} = nx Using Progression pr;
  Transfinite Line{8+fshift} = ny Using Progression pr;  
  Transfinite Line{9+fshift} = ny Using Progression pr;    
  Transfinite Line{10+fshift} = ny Using Progression pr;
  Transfinite Line{11+fshift} = ny Using Progression pr;
  Transfinite Line{12+fshift} = ny Using Progression pr;
  Transfinite Line{13+fshift} = ny Using Progression pr;
  Transfinite Line{14+fshift} = ny Using Progression pr;


// Definição da superfície 

  Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
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

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1,2,3,4,5,6,7,8,9,10,11,12};
  Physical Line("top") = {14,15,16,17,18,19,20,21,22,23,24};
  Physical Line("right") = {25};
  Physical Line("left") = {13};
  
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
  Physical Point("P31") = {31};
  Physical Point("P32") = {32};

  Physical Line("f14") = {14+fshift};      
  Physical Point("P33") = {33};
  Physical Point("P34") = {34};
  
  Coherence Mesh;


