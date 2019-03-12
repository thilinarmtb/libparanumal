r = DefineNumber[0.1];
Point(1) = {-1, -1, 0, r};
Point(2) = { 1, -1, 0, r};
Point(3) = {1, 1, 0, r};
Point(4) = {-1, 1, 0, r};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1,2,3,4};
Plane Surface(6) = {6};
Recombine Surface {6};
Physical Line("Inflow",2) = {1,2,3,4};
Physical Surface("Domain",9) = {6};