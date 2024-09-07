Include "../build/Data4Geo.txt";

Num_L = num_elem_L;
Prog_L =1;

Num_H = num_elem_H;
Prog_H =1;

Point(1) = {0, H ,0,1};
Point(2) = {0, -H,0,1};
Point(3) = {L, -H,0,1};
Point(4) = {L,  H,0,1};

Line(1) = {1,2};   Transfinite Curve {1} = Num_H Using Progression Prog_H;
Line(2) = {2,3};   Transfinite Curve {2} = Num_L Using Progression Prog_L;
Line(3) = {3,4};   Transfinite Curve {3} = Num_H Using Progression Prog_H;
Line(4) = {4,1};   Transfinite Curve {4} = Num_L Using Progression Prog_L;

Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

//PHYSICAL QUANTITIES

Physical Surface(1) = {1};
Physical Curve(1) = {1};   //left
Physical Curve(2) = {3};    //right