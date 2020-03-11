
// Gmsh project created on Fri Mar 09 15:42:42 2018
SetFactory("OpenCASCADE");
//R�alisation d'un bureau/table rectangulaire avec 4 pieds rectagulaires
//ltable=1;
//Ltable=0.5;
//etable=
//epied=50;
//longpied=40
//hpied=800;

lc=10.0e-3;

Point(1) = {0, 0, 0, lc};
Point(2) = {0, 40e-3, 0, lc};
Point(3) = {50e-3, 0, 0, lc};
Point(4) = {50e-3, 40e-3, 0, lc};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 2};
//+
Line(4) = {2, 1};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Surface(1) = {1};
//+
Extrude {0, 0, 600e-3} {
  Surface{1}; 
}


//+
Point(9) = {950e-3, 0, 0, lc};
//+
Point(10) = {1000e-3, 0, 0, lc};
//+
Point(11) = {950e-3, 40e-3, 0, lc};
//+
Point(12) = {1000e-3, 40e-3, 0, lc};
//+
Line(13) = {9, 10};
//+
Line(14) = {10, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {11, 9};
//+
Line Loop(8) = {16, 13, 14, 15};
//+
Surface(7) = {8};
//+
Extrude {0, 0, 600e-3} {
  Surface{7}; 
}
//+
Point(17) = {0, 460e-3, 0, lc};
//+
Point(18) = {50e-3, 460e-3, 0, lc};
//+
Point(19) = {50e-3, 500e-3, 0, lc};
//+
Point(20) = {950e-3, 460e-3, 0, lc};
//+
Point(21) = {1000e-3, 460e-3, 0, lc};
//+
Point(22) = {1000e-3, 500e-3, 0, lc};
//+
Point(23) = {950e-3, 500e-3, 0, lc};
//+
Point(24) = {0, 500e-3, 0, lc};
//+
Line(25) = {17, 18};
//+
Line(26) = {18, 19};
//+
Line(27) = {19, 24};
//+
Line(28) = {24, 17};
//+
Line(29) = {20, 21};
//+
Line(30) = {21, 22};
//+
Line(31) = {22, 23};
//+
Line(32) = {23, 20};
//+
Line Loop(15) = {25, 26, 27, 28};
//+
Surface(13) = {15};
//+
Line Loop(17) = {30, 31, 32, 29};
//+
Surface(14) = {17};
//+
Extrude {0, 0, 600e-3} {
  Surface{13}; Surface{14}; 
}
//+
Line(57) = {5, 15};
//+
Line(58) = {15, 37};
//+
Line(59) = {37, 35};
//+
Line(60) = {35, 5};
//+
Line Loop(29) = {59, 60, 57, 58};
//+
Surface(25) = {29};
//+
Extrude {0, 0, 20e-3} {
  Surface{25}; 
}

