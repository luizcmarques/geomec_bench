
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

  Point(1+dshift) = {26.70068359375, 10, 0, 1e+22};
  Point(2+dshift) = {29.08984375, 6.23486328125, 0, 1e+22};
  Point(3+dshift) = {93.11572265625, 10, 0, 1e+22}; 
  Point(4+dshift) = {94.50732421875, 0, 0, 1e+22}; 
  Point(5+dshift) = {122.5947265625, 10, 0, 1e+22}; 
  Point(6+dshift) = {127.7412109375, 0, 0, 1e+22}; 
  Point(7+dshift) = {193.3662109375, 9.9541015625, 0, 1e+22}; 
  Point(8+dshift) = {200, 1.23388671875, 0, 1e+22}; 
  Point(9+dshift) = {154.40966796875, 10, 0, 1e+22}; 
  Point(10+dshift) = {160.11865234375, 0, 0, 1e+22}; 
  Point(11+dshift) = {56.56591796875, 10, 0, 1e+22}; 
  Point(12+dshift) = {60.287109375, 0, 0, 1e+22}; 
  Point(13+dshift) = {120.54345703125, 10, 0, 1e+22}; 
  Point(14+dshift) = {125.86328125,  0, 0, 1e+22}; 
  Point(15+dshift) = {190.78125, 4.248046875, 0, 1e+22}; 
  Point(16+dshift) = {193.1943359375, 0, 0, 1e+22}; 
  Point(17+dshift) = {51.29345703125, 2.9326171875, 0, 1e+22}; 
  Point(18+dshift) = {53.435546875, 0, 0, 1e+22}; 

  Point(103) = {120.7436378,9.623708010335909, 0, 1e+22}; //Intersection

  Point(19+dshift) = {120.6572265625, 10, 0, 1e+22}; 
  Point(20+dshift) = {122.95361328125, 0, 0, 1e+22}; 
  Point(21+dshift) = {138.86767578125, 6.60302734375, 0, 1e+22}; 
  Point(22+dshift) = {141.1396484375, 0, 0, 1e+22}; 
  Point(23+dshift) = {68.560546875, 6.646484375, 0, 1e+22}; 
  Point(24+dshift) = {72.30810546875, 0, 0, 1e+22}; 
  Point(25+dshift) = {7.595703125, 10, 0, 1e+22}; 
  Point(26+dshift) = {11.11376953125, 0, 0, 1e+22}; 
  Point(27+dshift) = {189.45458984375, 10, 0, 1e+22}; 
  Point(28+dshift) = {193.30615234375, 0, 0, 1e+22}; 
  Point(29+dshift) = {197.37255859375, 10, 0, 1e+22}; 
  Point(30+dshift) = {199.87548828125, 0, 0, 1e+22}; 

  Point(104) = {199.35416716586613,2.08284362915751, 0, 1e+22}; //Intersection

  Point(31+dshift) = {159.3076171875, 10, 0, 1e+22}; 
  Point(32+dshift) = {159.8720703125, 6.67578125, 0, 1e+22}; 
  Point(33+dshift) = {66.70068359375, 10, 0, 1e+22}; 
  Point(34+dshift) = {67.61865234375,  0, 0, 1e+22}; 
  Point(35+dshift) = {117.12451171875, 10, 0, 1e+22}; 
  Point(36+dshift) = {119.69580078125, 0, 0, 1e+22}; 
  Point(37+dshift) = {150.11669921875, 10, 0, 1e+22}; 
  Point(38+dshift) = {155.00732421875, 0, 0, 1e+22}; 
  Point(39+dshift) = {130.556640625, 10, 0, 1e+22}; 
  Point(40+dshift) = {133.2509765625,  0, 0, 1e+22}; 
  Point(41+dshift) = {70.60595703125, 10, 0, 1e+22}; 
  Point(42+dshift) = {72.8486328125, 0, 0, 1e+22}; 
  Point(43+dshift) = {52.7529296875, 10, 0, 1e+22}; 
  Point(44+dshift) = {54.4169921875, 0, 0, 1e+22}; 
  Point(45+dshift) = {0.5771484375, 10, 0, 1e+22}; 
  Point(46+dshift) = {4.5419921875, 0, 0, 1e+22}; 
  
  Point(47+dshift) = {4.32177734375, 10, 0, 1e+22}; 
  Point(48+dshift) = {9.75439453125, 0, 0, 1e+22}; 
  Point(49+dshift) = {5.6865234375, 10, 0, 1e+22}; 
  Point(50+dshift) = {8.70703125, 0, 0, 1e+22}; 
  Point(100) = {7.395494962993421, 4.342105263157896, 0, 1e+22}; //Intersection

  Point(51+dshift) = {63.36669921875, 10, 0, 1e+22}; 
  Point(52+dshift) = {68.806640625, 0, 0, 1e+22}; 
  Point(101) = {67.37748854191973, 2.627146097, 0, 1e+22}; //Intersection

  Point(53+dshift) = {147.62890625, 6.15087890625, 0, 1e+22}; 
  Point(54+dshift) = {148.50830078125, 0, 0, 1e+22}; 
  Point(55+dshift) = {76.34716796875, 9.419921875, 0, 1e+22}; 
  Point(56+dshift) = {79.30419921875, 0, 0, 1e+22}; 
  Point(57+dshift) = {119.36376953125, 7.1572265625, 0, 1e+22}; 
  Point(58+dshift) = {123.23486328125, 0, 0, 1e+22};

  Point(102) = {122.74609271220551,0.9036830223, 0, 1e+22}; //Intersection

  Point(59+dshift) = {139.31201171875, 0, 0, 1e+22}; 
  Point(60+dshift) = {139.408203125, 4.21630859375, 0, 1e+22}; 

// Fronteiras - Domínio Omega  

  Line(1) = {1,50}; 
  Line(2) = {50,54};
  Line(3) = {54,52};
  Line(4) = {52,30};
  Line(5) = {30,22};
  Line(6) = {22,48};
  Line(7) = {48,16};
  Line(8) = {16,38};
  Line(9) = {38,56};
  Line(10) = {56,28};
  Line(11) = {28,46};
  Line(12) = {46,60};
  Line(13) = {60,8};
  Line(14) = {8,40};
  Line(15) = {40,24};
  Line(16) = {24,62};
  Line(17) = {62,18};
  Line(18) = {18,10};
  Line(19) = {10,44};
  Line(20) = {44,63};
  Line(21) = {63,26};
  Line(22) = {26,58};
  Line(23) = {58,42};
  Line(24) = {42,14};
  Line(25) = {14,20};
  Line(26) = {32,34};
  Line(27) = {20,32};
  Line(28) = {34,2};

  Line(29) = {2,12};
  Line(30) = {12,3};
 
  Line(31) = {3,33};
  Line(32) = {33, 31};
  Line(33) = {31, 35};
  Line(34) = {35, 13};
  Line(35) = {13, 41};
  Line(36) = {41, 43};
  Line(37) = {43, 9};
  Line(38) = {9, 23}; 
  Line(39) = {23, 17};
  Line(40) = {17, 39};

  Line(41) = {39, 7}; 
  Line(42) = {7, 45};
  Line(43) = {45, 37};
  Line(44) = {37, 55};
  Line(45) = {55, 15};
  Line(46) = {15, 47};
  Line(47) = {47, 5};
  Line(48) = {5, 29};
  Line(49) = {29, 53};
  Line(50) = {53, 51};
  Line(51) = {51, 49};
  Line(52) = {49, 4};
  
  Line(53) = {4,1};

//Fratura

  fshift = 53;

  Line(1+fshift) = {1+dshift,2+dshift};
  Line(2+fshift) = {3+dshift,4+dshift};
  Line(3+fshift) = {5+dshift,6+dshift};
  
  Line(4+fshift+100) = {7+dshift,104};
  Line(4+fshift+101) = {104,8+dshift};

  Line(5+fshift) = {9+dshift,10+dshift};
  Line(6+fshift) = {11+dshift,12+dshift};
  
  Line(7+fshift+102) = {13+dshift,103};
  Line(7+fshift+103) = {103,14+dshift};
  
  
  Line(8+fshift) = {15+dshift,16+dshift};
  Line(9+fshift) = {17+dshift,18+dshift};

  Line(10+fshift+104) = {19+dshift,103};
  Line(10+fshift+105) = {102,103};
  Line(10+fshift+106) = {102,20+dshift};

  Line(11+fshift) = {21+dshift,22+dshift};
  Line(12+fshift) = {23+dshift,24+dshift};
  Line(13+fshift) = {25+dshift,26+dshift};
  Line(14+fshift) = {27+dshift,28+dshift};
  
  Line(15+fshift+107) = {29+dshift,104};
  Line(15+fshift+108) = {104,30+dshift};

  Line(16+fshift) = {31+dshift,32+dshift};
  
  Line(17+fshift+109) = {33+dshift,101};
  Line(17+fshift+110) = {101,34+dshift};

  Line(18+fshift) = {35+dshift,36+dshift};
  Line(19+fshift) = {37+dshift,38+dshift};
  Line(20+fshift) = {39+dshift,40+dshift};
  Line(21+fshift) = {41+dshift,42+dshift};
  Line(22+fshift) = {43+dshift,44+dshift};
  Line(23+fshift) = {45+dshift,46+dshift};
  
  Line(24+fshift+111) = {47+dshift,100};
  Line(24+fshift+112) = {100,48+dshift};

  Line(25+fshift+113) = {49+dshift,100};
  Line(25+fshift+114) = {100,50+dshift};
  
  Line(26+fshift+115) = {51+dshift,101};
  Line(26+fshift+116) = {101,52+dshift};

  Line(27+fshift) = {53+dshift,54+dshift};
  Line(28+fshift) = {55+dshift,56+dshift};
  
  Line(29+fshift+117) = {57+dshift,102};
  Line(29+fshift+118) = {102,58+dshift};

  Line(30+fshift) = {59+dshift,60+dshift};


  Transfinite Line{5,13,14,22,24,25,33,36,41,42,47,48,53} = nx Using Progression pr;
  Transfinite Line{1,2,3,4,6,7,8,9,10,11,12,15,16,17,18,19,20,21,23,26,29,30,31,32,34,35,37,38,39,40,43,44,45,46,49,50,51,52} = ny Using Progression pr;

  Transfinite Line{1+fshift} = ny Using Progression pr;
  Transfinite Line{2+fshift} = ny Using Progression pr;
  Transfinite Line{3+fshift} = ny Using Progression pr;    
  Transfinite Line{4+fshift+100} = ny Using Progression pr;
  Transfinite Line{4+fshift+101} = ny Using Progression pr;
  Transfinite Line{5+fshift} = ny Using Progression pr;
  Transfinite Line{6+fshift} = ny Using Progression pr;
  Transfinite Line{7+fshift+102} = nx Using Progression pr;
  Transfinite Line{7+fshift+103} = nx Using Progression pr;
  Transfinite Line{8+fshift} = ny Using Progression pr;  
  Transfinite Line{9+fshift} = ny Using Progression pr;    
  Transfinite Line{10+fshift+104} = ny Using Progression pr;
  Transfinite Line{10+fshift+105} = ny Using Progression pr;
  Transfinite Line{10+fshift+106} = ny Using Progression pr;
  Transfinite Line{11+fshift} = ny Using Progression pr;
  Transfinite Line{12+fshift} = ny Using Progression pr;
  Transfinite Line{13+fshift} = ny Using Progression pr;
  Transfinite Line{14+fshift} = ny Using Progression pr;
  Transfinite Line{15+fshift+107} = ny Using Progression pr;
  Transfinite Line{15+fshift+108} = ny Using Progression pr;
  Transfinite Line{16+fshift} = ny Using Progression pr;
  Transfinite Line{17+fshift+109} = ny Using Progression pr;
  Transfinite Line{17+fshift+110} = ny Using Progression pr;    
  Transfinite Line{18+fshift} = ny Using Progression pr;
  Transfinite Line{19+fshift} = ny Using Progression pr;
  Transfinite Line{20+fshift} = ny Using Progression pr;
  Transfinite Line{21+fshift} = nx Using Progression pr;
  Transfinite Line{22+fshift} = ny Using Progression pr;  
  Transfinite Line{23+fshift} = ny Using Progression pr;    
  Transfinite Line{24+fshift+111} = ny Using Progression pr;
  Transfinite Line{24+fshift+112} = ny Using Progression pr;
  Transfinite Line{25+fshift+113} = ny Using Progression pr;
  Transfinite Line{25+fshift+114} = ny Using Progression pr;
  Transfinite Line{26+fshift+115} = ny Using Progression pr;
  Transfinite Line{26+fshift+116} = ny Using Progression pr;
  Transfinite Line{27+fshift} = ny Using Progression pr;
  Transfinite Line{28+fshift} = ny Using Progression pr;
  Transfinite Line{29+fshift+117} = ny Using Progression pr;
  Transfinite Line{29+fshift+118} = ny Using Progression pr;
  Transfinite Line{30+fshift} = ny Using Progression pr;


// Definição da superfície 

  Line Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53};
  Plane Surface(1) = {1};

  Line{1+fshift} In Surface{1};
  Line{2+fshift} In Surface{1};
  Line{3+fshift} In Surface{1};
  Line{4+fshift+100} In Surface{1};
  Line{4+fshift+101} In Surface{1};
  Line{5+fshift} In Surface{1};
  Line{6+fshift} In Surface{1};
  Line{7+fshift+102} In Surface{1};
  Line{7+fshift+103} In Surface{1};
  Line{8+fshift} In Surface{1};
  Line{9+fshift} In Surface{1};
  Line{10+fshift+104} In Surface{1};
  Line{10+fshift+105} In Surface{1};
  Line{10+fshift+106} In Surface{1};
  Line{11+fshift} In Surface{1};
  Line{12+fshift} In Surface{1};
  Line{13+fshift} In Surface{1};
  Line{14+fshift} In Surface{1};
  Line{15+fshift+107} In Surface{1};
  Line{15+fshift+108} In Surface{1};
  Line{16+fshift} In Surface{1};
  Line{17+fshift+109} In Surface{1};
  Line{17+fshift+110} In Surface{1};
  Line{18+fshift} In Surface{1};
  Line{19+fshift} In Surface{1};
  Line{20+fshift} In Surface{1};
  Line{21+fshift} In Surface{1};
  Line{22+fshift} In Surface{1};
  Line{23+fshift} In Surface{1};
  Line{24+fshift+111} In Surface{1};
  Line{24+fshift+112} In Surface{1};
  Line{25+fshift+113} In Surface{1};
  Line{25+fshift+114} In Surface{1};
  Line{26+fshift+115} In Surface{1};
  Line{26+fshift+116} In Surface{1};
  Line{27+fshift} In Surface{1};
  Line{28+fshift} In Surface{1};
  Line{29+fshift+117} In Surface{1};
  Line{29+fshift+118} In Surface{1};
  Line{30+fshift} In Surface{1};

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};
  Physical Line("top") = {31,32,33,34,35,36,37,38,40,41,42,43,44,45,46,47,48,49,50,51,52,53};
  Physical Line("right") = {29,30};
  Physical Line("left") = {54};
  
  Physical Line("f1") = {1+fshift};
  Physical Point("P5") = {5};
  Physical Point("P6") = {6};

  Physical Line("f2") = {2+fshift};
  Physical Point("P7") = {7};
  Physical Point("P8") = {8};

  Physical Line("f3") = {3+fshift};
  Physical Point("P9") = {9};
  Physical Point("P10") = {10};

  Physical Line("f4") = {4+fshift+100,4+fshift+101};
  Physical Point("P11") = {11};
  Physical Point("P12") = {12};

  Physical Line("f5") = {5+fshift};
  Physical Point("P13") = {13};
  Physical Point("P14") = {14};

  Physical Line("f6") = {6+fshift};
  Physical Point("P15") = {15};
  Physical Point("P16") = {16};

  Physical Line("f7") = {7+fshift+102,7+fshift+103};
  Physical Point("P17") = {17};
  Physical Point("P18") = {18};

  Physical Line("f8") = {8+fshift};
  Physical Point("P19") = {19};
  Physical Point("P20") = {20};

  Physical Line("f9") = {9+fshift};
  Physical Point("P21") = {21};
  Physical Point("P22") = {22};

  Physical Line("f10") = {10+fshift+104,10+fshift+105,10+fshift+106};
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

  Physical Line("f15") = {15+fshift+107,15+fshift+108};      
  Physical Point("P33") = {33};
  Physical Point("P34") = {34};

  Physical Line("f16") = {16+fshift};
  Physical Point("P35") = {35};
  Physical Point("P36") = {36};

  Physical Line("f17") = {17+fshift+109,17+fshift+110};
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

  Physical Line("f21") = {21+fshift};
  Physical Point("P45") = {45};
  Physical Point("P46") = {46};

  Physical Line("f22") = {22+fshift};      
  Physical Point("P47") = {47};
  Physical Point("P48") = {48};

  Physical Line("f23") = {23+fshift};      
  Physical Point("P49") = {49};
  Physical Point("P50") = {50};

  Physical Line("f24") = {24+fshift+111,24+fshift+112};      
  Physical Point("P51") = {51};
  Physical Point("P52") = {52};

  Physical Line("f25") = {25+fshift+113,25+fshift+114};      
  Physical Point("P53") = {53};
  Physical Point("P54") = {54};

  Physical Line("f26") = {26+fshift+115,26+fshift+116};      
  Physical Point("P55") = {55};
  Physical Point("P56") = {56};

  Physical Line("f27") = {27+fshift};      
  Physical Point("P57") = {57};
  Physical Point("P58") = {58};  

  Physical Line("f28") = {28+fshift};      
  Physical Point("P59") = {59};
  Physical Point("P60") = {60}; 

  Physical Line("f29") = {29+fshift+117,29+fshift+118};      
  Physical Point("P61") = {61};
  Physical Point("P62") = {62};   

  Physical Line("f30") = {30+fshift};      
  Physical Point("P63") = {63};
  Physical Point("P64") = {64};  
  Coherence Mesh;


