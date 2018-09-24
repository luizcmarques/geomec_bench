//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-0.8, -0.4, 0.2, 1.5, 1.1, 0};
//+
Characteristic Length {1, 2, 3, 4} = 2;
//+

Transfinite Line {1,2,3,4} = 1 Using Progression 1;

Transfinite Surface {1};

Recombine Surface {1};

//+
Physical Surface("domain") = {1};
