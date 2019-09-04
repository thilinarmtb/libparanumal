 res = DefineNumber[0.1];
 
 xn = DefineNumber[100];
 yn = DefineNumber[ 10];

 xmax = DefineNumber[10.0]; 
 xmin = DefineNumber[ 0.0]; 
 ymax = DefineNumber[ 1.0]; 
 ymin = DefineNumber[ 0.0]; 

 Point(1) = {xmin, ymin, 0, res};
 Point(2) = {xmax, ymin, 0, res};
 Point(3) = {xmax, ymax, 0, res};
 Point(4) = {xmin, ymax, 0, res};
 
 Line(1) = {1, 2};
 Line(2) = {2, 3};
 Line(3) = {3, 4};
 Line(4) = {4, 1};
 
 Transfinite Line {1, 3} = (xn+1) Using Progression 1;
 Transfinite Line {2, 4} = (yn+1) Using Progression 1;

 Line Loop(9) = {1, 2, 3, 4};
 Plane Surface(9) = {9};
 Transfinite Surface {9} = {1,2,3,4};
 
 Physical Surface("Domain",9) = {9};
 Physical Line("Inflow",2)    = {4};
 Physical Line("Wall",1)      = {1,3};
 Physical Line("Outflow",3)   = {2};
 Recombine Surface {9};

