// Gmsh project created on Tue Feb  7 12:11:26 2023
//+
lc1=1;
//+
lc2=1;
//+
lc3=1;
//+
Point(1) = {0, 0.0, 0.0, lc3};
//+
Point(2) = {200, 0.0, 0.0, lc3};
//+
Point(3) = {200, 150, 0.0, lc1};
//+
Point(4) = {0.0, 150, 0.0, lc1};
//+
Point(5) = {200, 114, 0.0, lc2};
//+
Point(6) = {0.0, 114, 0.0, lc2};
//+
Point(7) = {0, 114, 0.0, lc2};
//+
Point(8) = {200, 127, 0.0, lc2};
//+
Point(9) = {0.0, 127, 0.0, lc2};
//+
Point(10) = {200, 147.5, 0.0, lc1};
//+
Point(11) = {0.0, 147.5, 0.0, lc1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 8};
//+
Line(4) = {8, 10};
//+
Line(5) = {10, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 11};
//+
Line(8) = {11, 9};
//+
Line(9) = {9, 6};
//+
Line(10) = {6, 1};
//+
Line(11) = {5, 6};
//+
Line(12) = {8, 9};
//+
Line(13) = {10, 11};
//+
Line Loop(1) = {1, 2, 11, 10};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {11, -9, -12, -3};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {12, -8, -13, -4};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {13, -7, -6, -5};
//+
Plane Surface(4) = {4};
//+
Physical Line(1) = {1};
//+
Physical Line(2) = {2, 3, 4, 5, 7, 8, 9, 10};
//+
Physical Surface(3) = {1};
//+
Physical Surface(4) = {2};
//+
Physical Surface(5) = {3};
//+
Physical Surface(6) = {4};
//+
Extrude {0, 0, 200} {
  Line{10}; Line{1}; Line{2}; Surface{1}; Surface{2}; Surface{3}; Surface{4}; Line{12}; Line{11}; Line{3}; Line{4}; Line{13}; Line{5}; Line{8}; Line{9}; Point{4}; Line{6}; Line{7};
}
//+
Physical Volume(7) = {1};
//+
Physical Volume(8) = {2};
//+
Physical Volume(9) = {3};
//+
Physical Volume(10) = {4};
