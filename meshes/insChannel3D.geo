s = 0.1;

xmin = DefineNumber[0.0]; 
xmax = DefineNumber[10.0]; 

ymin = DefineNumber[0.0]; 
ymax = DefineNumber[0.5]; 

zmin = DefineNumber[-0.5]; 
zmax = DefineNumber[ 0.5]; 

Point(1) = {xmin,ymax, zmin, s};
Point(2) = {xmin,ymax, zmax, s};
Point(3) = {xmin,ymin, zmax, s};
Point(4) = {xmin,ymin, zmin, s};

Point(5) = {xmax,ymax, zmin, s};
Point(6) = {xmax,ymax, zmax, s};
Point(7) = {xmax,ymin, zmax, s};
Point(8) = {xmax,ymin, zmin, s};

Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};
Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 5};
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Line Loop(3) = {9, 5, -10, -1};
Plane Surface(3) = {3};
Line Loop(4) = {10, 6, -11, -2};
Plane Surface(4) = {4};
Line Loop(5) = {11, 7, -12, -3};
Plane Surface(5) = {5};
Line Loop(6) = {9, -8, -12, 4};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Physical Surface("Wall", 1) = {3,4,5,6};
Physical Surface("Wall", 1) = {3};
Physical Surface("Yslip", 5) = {5};
Physical Surface("Zslip", 6) = {4,6};

Physical Surface("Inflow", 2) = {1};
Physical Surface("Outflow", 3) = {2};
Physical Volume("Domain", 9) = {1};
