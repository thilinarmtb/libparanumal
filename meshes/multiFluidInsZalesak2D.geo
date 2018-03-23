
r = DefineNumber[0.025];

Point(1) = {0, 0, 0, r};
Point(2) = {1, 0, 0, r};
Point(3) = {1, 1, 0, r};
Point(4) = {0, 1, 0, r};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Physical Surface("Interior",9) = {6};
Physical Line("Inflow",2) = {1, 2, 3, 4};
