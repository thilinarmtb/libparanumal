r0 = DefineNumber[0.25];
r1 = DefineNumber[0.20];
r2 = DefineNumber[0.10];
r3 = DefineNumber[0.05];

// Domain Boundaries
xmin  = DefineNumber[00.0];
xmax  = DefineNumber[32.0];
ymax  = DefineNumber[16.0];
ymin  = DefineNumber[00.0];

// Cylinder coordinates
xcy0  = DefineNumber[08.0];
ycy0  = DefineNumber[08.0];
R     = DefineNumber[00.5];
xcy1  = xcy0 + R*Cos(45*Pi/180);
xcy2  = xcy0 - R*Cos(45*Pi/180);
xcy3  = xcy0 - R*Cos(45*Pi/180);
xcy4  = xcy0 + R*Cos(45*Pi/180);
//
ycy1  = ycy0 + R*Sin(45*Pi/180);
ycy2  = ycy0 + R*Sin(45*Pi/180);
ycy3  = ycy0 - R*Sin(45*Pi/180);
ycy4  = ycy0 - R*Sin(45*Pi/180);
// A box around cylinder
xbmax  = DefineNumber[10.0];
xbmin  = DefineNumber[06.0];
ybmax  = DefineNumber[10.0];
ybmin  = DefineNumber[06.0];




// Domain
Point(1) = {xmin, ymin, 0, r0};
Point(2) = {xmax, ymin, 0, r0};
Point(3) = {xmax, ymax, 0, r0};
Point(4) = {xmin, ymax, 0, r0};
// Cylinders
Point(5) = {xcy0, ycy0, 0, r3};
Point(6) = {xcy1, ycy1, 0, r3};
Point(7) = {xcy2, ycy2, 0, r3};
Point(8) = {xcy3, ycy3, 0, r3};
Point(9) = {xcy4, ycy4, 0, r3};
// Box
Point(10)={xbmin, ybmin, 0, r2};
Point(11)={xbmax, ybmin, 0, r2};
Point(12)={xbmax, ybmax, 0, r2};
Point(13)={xbmin, ybmax, 0, r2};

// Box Back
Point(14)={xmax, ybmin-2, 0, r1};
Point(15)={xmax, ybmax+2, 0, r1};

Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};

Line(5) = {10, 11};
Line(6) = {11, 12};
Line(7) = {12, 13};
Line(8) = {13, 10};
Line(9) = {12, 15};
Line(10) = {15, 14};
Line(11) = {14, 11};
Line(12) = {1, 2};
Line(13) = {2, 14};
Line(14) = {15, 3};
Line(15) = {3, 4};
Line(16) = {4, 1};
Line Loop(17) = {1, 2, 3, 4};
Line Loop(18) = {6, 7, 8, 5};
Plane Surface(19) = {17, 18};
Line Loop(20) = {6, 9, 10, 11};
Plane Surface(21) = {20};
Line Loop(22) = {12, 13, 11, -5, -8, -7, 9, 14, 15, 16};
Plane Surface(23) = {22};
Recombine Surface {23, 19, 21};

Coherence; 

Physical Line("Wall",1) = {1, 4, 3, 2};
Physical Line("Inflow",2) = {12, 16, 15};
Physical Line("Outflow",3) = {14, 10, 13};
Physical Surface("Domain",9) = {23, 21, 19};
