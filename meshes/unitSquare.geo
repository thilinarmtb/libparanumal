LC = 1.0;

Point(1) = {-1.0, -1.0, 0.0, LC};
Point(2) = { 1.0, -1.0, 0.0, LC};
Point(3) = { 1.0,  1.0, 0.0, LC};
Point(4) = {-1.0,  1.0, 0.0, LC};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Line{1} = 3;
Transfinite Line{2} = 3;
Transfinite Line{3} = 3;
Transfinite Line{4} = 3;

Recombine Surface{1};

RefineMesh;

Physical Line("Boundary") = {1, 2, 3, 4};
Physical Surface("Domain") = {1};
