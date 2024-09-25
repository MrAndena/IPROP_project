// ---- PARAMETERS ----

// RADIUS - DISTANCES - BOXES LENGTH:

r_emi= 30e-6; //[m]   radius of the emitter
r_col= 1e-3; // [m]   radius of the collector 

dist_emi_col = 0.02;    //  [m] distance emitter collector (starting from the perimiter)
dist_col_up_down = 0.01; //  [m] distance between the surface of the collector and the upper/bottom part of the domain 
dist_col_outlet = 0.015; //  [m] distance between the surface of the collector and the outlet 
dist_emi_inlet = 0.015; //   [m] distance between the surface of the emitter and the inlet


L_col_in = 2*r_col;       // [m] half the length of the inner square that wraps the collector (proportion to the radius)
L_col_out = 3.5*r_col;    // [m] half the length of the outer square that wraps the collector 

L_emi_in = 30*r_emi;       // [m] half the length of the inner square that wraps the emitter
L_emi_out = 60*r_emi;      // [m] half the length of the outer square that wraps the emitter
                          // NB: for what concern the emitter, there is asecond other box that wraps it; it has the dimension of L_col_out

// COLLECTOR TRANSFINITE VALUES:

Num_col_circ = 26;      //number of points of the transfinite lines on the collector circunference
Prog_col_circ = 1;      //progression number on collector circunference

Num_col_in_box = Num_col_circ;    //number of points of the transfinite lines on the collector inner box
Prog_col_in_box = 1;              //progression number on collector inner box
Num_col_in_L = 31;                //number of points of the transfinite lines on the collector inner box obliquos edge
Prog_col_in_L = 1;                //progression number on collector inner box obl. edge

// outer box that contains emitter and collector transfinite value
Num_box = Num_col_circ;
Prog_box = 1;
Num_box_L = 26;
Prog_box_L =1;

// EMITTER TRANSFINITE VALUES:

Num_emi_circ = Num_box;      //number of points of the transfinite lines on the emitter circunference
Prog_emi_circ = 1;           //progression number on the emitter circunference

Num_emi_in_box = Num_box;    //number of points of the transfinite lines on the emitter inner box
Prog_emi_in_box = 1;         //progression number on emitter inner box
Num_emi_in_L = 26;           //number of points of the transfinite lines on the inner emitter box obl. edges
Prog_emi_in_L = 1;           //progression number on emitter inner box obl. edge

Num_emi_out_box = Num_box;   //number of points of the transfinite lines on the emitter outer box
Prog_emi_out_box = 1;        //progression number on emitter outer box
Num_emi_out_L = 26;          //number of points of the transfinite lines on the outer emitter box obl. edges
Prog_emi_out_L = 1;          //progression number on emitter outer box obl. edge

// recall that the emitter has another outher box, see in section of collector

// DOMAIN TRANSFINITE VALUES: 

Num_dom_left_left =26;      //number of points of the transfinite lines on the left edge of the domain
Prog_dom_left_left =1;      //progression number on the left edge of the domain

Num_dom_left_up =26;        //number of points of the transfinite lines on the left up corner of the domain
Prog_dom_left_up=1;         //progression number on the left up corner of the domain

Num_dom_right_up =26;       //number of points of the transfinite lines on the right up corner of the domain
Prog_dom_right_up =1;       //progression number on the right up corner of the domain

Num_dom_center =26;         //number of points of the transfinite lines on the central edge of the domain
Prog_dom_center=1;          //progression number on the central edge of the domain

// parameters for the progression of the transfinite lines, I need two since i could have different orientation in the SURFACES
A = 1.00;
B = 1.00;

// maybe can be useful
mesh_ref_1 = 1.0;

// ---- POINTS ----

// COLLECTOR AND INNER BOX COLLECTOR

Point(1) = {r_col , 0, 0, mesh_ref_1 };    // center collector
Point(2) = {r_col + r_col/Sqrt(2) , r_col/Sqrt(2), 0, mesh_ref_1 };
Point(3) = {r_col - r_col/Sqrt(2) , r_col/Sqrt(2), 0, mesh_ref_1 };
Point(4) = {r_col - r_col/Sqrt(2) , -r_col/Sqrt(2), 0, mesh_ref_1};
Point(5) = {r_col + r_col/Sqrt(2) , -r_col/Sqrt(2), 0, mesh_ref_1};

Point(6) = {r_col + L_col_in, L_col_in, 0, mesh_ref_1};
Point(7) = {r_col - L_col_in, L_col_in, 0, mesh_ref_1};
Point(8) = {r_col - L_col_in, -L_col_in, 0, mesh_ref_1};
Point(9) = {r_col + L_col_in, -L_col_in, 0, mesh_ref_1};

// EMITTER AND INNER BOX EMITTER 

Point(10) = {-dist_emi_col - r_emi, 0, 0, mesh_ref_1}; // center emitter
Point(11) = {-dist_emi_col - r_emi + r_emi/Sqrt(2), r_emi/Sqrt(2), 0, mesh_ref_1};
Point(12) = {-dist_emi_col - r_emi - r_emi/Sqrt(2), r_emi/Sqrt(2), 0, mesh_ref_1};
Point(13) = {-dist_emi_col - r_emi - r_emi/Sqrt(2), -r_emi/Sqrt(2), 0, mesh_ref_1};
Point(14) = {-dist_emi_col - r_emi + r_emi/Sqrt(2), -r_emi/Sqrt(2), 0, mesh_ref_1};

Point(15) = {-dist_emi_col - r_emi + L_emi_in, L_emi_in, 0, mesh_ref_1};
Point(16) = {-dist_emi_col - r_emi - L_emi_in, L_emi_in, 0, mesh_ref_1};
Point(17) = {-dist_emi_col - r_emi - L_emi_in, -L_emi_in, 0, mesh_ref_1};
Point(18) = {-dist_emi_col - r_emi + L_emi_in, -L_emi_in, 0, mesh_ref_1};

// MAIN DOMAIN 

Point(19) = {2*r_col + dist_col_outlet, r_col + dist_col_up_down, 0, mesh_ref_1};
Point(20) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, r_col + dist_col_up_down, 0, mesh_ref_1};
Point(21) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, -r_col - dist_col_up_down, 0, mesh_ref_1};
Point(22) = {2*r_col + dist_col_outlet, -r_col - dist_col_up_down, 0, mesh_ref_1};

// OUTER BOXES AND MORE DOMAIN POINTS
// coll
Point(23) = {r_col + L_col_out, L_col_out, 0, mesh_ref_1};
Point(24) = {r_col - L_col_out, L_col_out, 0, mesh_ref_1};
Point(25) = {r_col - L_col_out, -L_col_out, 0, mesh_ref_1};
Point(26) = {r_col + L_col_out, -L_col_out, 0, mesh_ref_1};
// emi 1
Point(27) = {-dist_emi_col - r_emi + L_emi_out, L_emi_out, 0, mesh_ref_1};
Point(28) = {-dist_emi_col - r_emi - L_emi_out, L_emi_out, 0, mesh_ref_1};
Point(29) = {-dist_emi_col - r_emi - L_emi_out, -L_emi_out, 0, mesh_ref_1};
Point(30) = {-dist_emi_col - r_emi + L_emi_out, -L_emi_out, 0, mesh_ref_1};
// emi 2
Point(31) = {-dist_emi_col - r_emi + L_col_out, L_col_out, 0, mesh_ref_1};
Point(32) = {-dist_emi_col - r_emi - L_col_out, L_col_out, 0, mesh_ref_1};
Point(33) = {-dist_emi_col - r_emi - L_col_out, -L_col_out, 0, mesh_ref_1}; 
Point(34) = {-dist_emi_col - r_emi + L_col_out, -L_col_out, 0, mesh_ref_1};

// POINTS ON THE PERIMITER OF THE RECTANGULAR DOMAIN

Point(35) = {r_col + L_col_out,r_col + dist_col_up_down, 0, mesh_ref_1};
Point(36) = {r_col - L_col_out,r_col + dist_col_up_down, 0, mesh_ref_1};
Point(37) = {-dist_emi_col - r_emi + L_col_out, r_col + dist_col_up_down, 0, mesh_ref_1};
Point(38) = {-dist_emi_col - r_emi - L_col_out, r_col + dist_col_up_down, 0, mesh_ref_1};
Point(39) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, L_col_out, 0, mesh_ref_1};
Point(40) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, -L_col_out, 0, mesh_ref_1}; 
Point(41) = {-dist_emi_col - r_emi - L_col_out,-r_col - dist_col_up_down, 0, mesh_ref_1};
Point(42) = {-dist_emi_col - r_emi + L_col_out,-r_col - dist_col_up_down, 0, mesh_ref_1};
Point(43) = {r_col - L_col_out,-r_col - dist_col_up_down, 0, mesh_ref_1};
Point(44) = {r_col + L_col_out,-r_col - dist_col_up_down, 0, mesh_ref_1};
Point(45) = {2*r_col + dist_col_outlet,-L_col_out, 0, mesh_ref_1};
Point(46) = {2*r_col + dist_col_outlet,L_col_out, 0, mesh_ref_1};
Point(47) = {-0.5*dist_emi_col, r_col + dist_col_up_down, 0, mesh_ref_1};
Point(48) = {-0.5*dist_emi_col, L_col_out, 0, mesh_ref_1};
Point(49) = {-0.5*dist_emi_col, -L_col_out, 0, mesh_ref_1};
Point(50) = {-0.5*dist_emi_col, -r_col - dist_col_up_down, 0, mesh_ref_1};


// ---- LINES ----

// COLLECTOR

Circle(1) = {2,1,3}; Transfinite Curve {1} = Num_col_circ Using Progression Prog_col_circ;
Circle(2) = {3,1,4}; Transfinite Curve {2} = Num_col_circ Using Progression Prog_col_circ;
Circle(3) = {4,1,5}; Transfinite Curve {3} = Num_col_circ Using Progression Prog_col_circ;
Circle(4) = {5,1,2}; Transfinite Curve {4} = Num_col_circ Using Progression Prog_col_circ;

Line(5) = {6,7};     Transfinite Curve {5} = Num_col_in_box Using Progression Prog_col_in_box;
Line(6) = {7,8};     Transfinite Curve {6} = Num_col_in_box Using Progression Prog_col_in_box;
Line(7) = {8,9};     Transfinite Curve {7} = Num_col_in_box Using Progression Prog_col_in_box;
Line(8) = {9,6};     Transfinite Curve {8} = Num_col_in_box Using Progression Prog_col_in_box;

Line(9) = {2,6};     Transfinite Curve {9} = Num_col_in_L Using Progression A*Prog_col_in_L;
Line(10) = {3,7};    Transfinite Curve {10} = Num_col_in_L Using Progression A*Prog_col_in_L;
Line(11) = {4,8};    Transfinite Curve {11} = Num_col_in_L Using Progression A*Prog_col_in_L;
Line(12) = {5,9};    Transfinite Curve {12} = Num_col_in_L Using Progression A*Prog_col_in_L;

// EMITTER

Circle(13) = {11,10,12};   Transfinite Curve {13} = Num_emi_circ Using Progression Prog_emi_circ;
Circle(14) = {12,10,13};   Transfinite Curve {14} = Num_emi_circ Using Progression Prog_emi_circ;
Circle(15) = {13,10,14};   Transfinite Curve {15} = Num_emi_circ Using Progression Prog_emi_circ;
Circle(16) = {14,10,11};   Transfinite Curve {16} = Num_emi_circ Using Progression Prog_emi_circ;

Line(17) = {15,16};        Transfinite Curve {17} = Num_emi_in_box Using Progression Prog_emi_in_box;
Line(18) = {16,17};        Transfinite Curve {18} = Num_emi_in_box Using Progression Prog_emi_in_box;
Line(19) = {17,18};        Transfinite Curve {19} = Num_emi_in_box Using Progression Prog_emi_in_box;
Line(20) = {18,15};        Transfinite Curve {20} = Num_emi_in_box Using Progression Prog_emi_in_box;

Line(21) = {11,15};        Transfinite Curve {21} = Num_emi_in_L Using Progression A*Prog_emi_in_L;
Line(22) = {12,16};        Transfinite Curve {22} = Num_emi_in_L Using Progression A*Prog_emi_in_L;
Line(23) = {13,17};        Transfinite Curve {23} = Num_emi_in_L Using Progression A*Prog_emi_in_L;
Line(24) = {14,18};        Transfinite Curve {24} = Num_emi_in_L Using Progression A*Prog_emi_in_L;

// DOMAIN

Line(25) = {23,24};    Transfinite Curve {25} = Num_box Using Progression Prog_box;   
Line(26) = {24,25};    Transfinite Curve {26} = Num_box Using Progression Prog_box; 
Line(27) = {25,26};    Transfinite Curve {27} = Num_box Using Progression Prog_box;
Line(28) = {26,23};    Transfinite Curve {28} = Num_box Using Progression Prog_box;
Line(29) = {6,23};    Transfinite Curve {29} = Num_box_L Using Progression Prog_box_L;
Line(30) = {7,24};    Transfinite Curve {30} = Num_box_L Using Progression Prog_box_L; 
Line(31) = {8,25};    Transfinite Curve {31} = Num_box_L Using Progression Prog_box_L; 
Line(32) = {9,26};    Transfinite Curve {32} = Num_box_L Using Progression Prog_box_L; 

Line(33) = {27,28};    Transfinite Curve {33} = Num_emi_out_box Using Progression Prog_emi_out_box;
Line(34) = {28,29};    Transfinite Curve {34} = Num_emi_out_box Using Progression Prog_emi_out_box; 
Line(35) = {29,30};    Transfinite Curve {35} = Num_emi_out_box Using Progression Prog_emi_out_box;
Line(36) = {30,27};    Transfinite Curve {36} = Num_emi_out_box Using Progression Prog_emi_out_box;
Line(37) = {15,27};    Transfinite Curve {37} = Num_emi_out_L Using Progression Prog_emi_out_L;
Line(38) = {16,28};    Transfinite Curve {38} = Num_emi_out_L Using Progression Prog_emi_out_L;
Line(39) = {17,29};    Transfinite Curve {39} = Num_emi_out_L Using Progression Prog_emi_out_L;
Line(40) = {18,30};    Transfinite Curve {40} = Num_emi_out_L Using Progression Prog_emi_out_L;
Line(41) = {31,32};    Transfinite Curve {41} = Num_box Using Progression Prog_box;
Line(42) = {32,33};    Transfinite Curve {42} = Num_box Using Progression Prog_box;
Line(43) = {33,34};    Transfinite Curve {43} = Num_box Using Progression Prog_box;
Line(44) = {34,31};    Transfinite Curve {44} = Num_box Using Progression Prog_box;

Line(45) = {27,31};    Transfinite Curve {45} = 26 Using Progression A*Prog_box_L;
Line(46) = {28,32};    Transfinite Curve {46} = 26 Using Progression A*Prog_box_L;
Line(47) = {29,33};    Transfinite Curve {47} = 26 Using Progression A*Prog_box_L;
Line(48) = {30,34};    Transfinite Curve {48} = 26 Using Progression A*Prog_box_L;

Line(51) = {44,22};    Transfinite Curve {51} = Num_dom_right_up Using Progression A*Prog_dom_right_up;
Line(52) = {26,45};    Transfinite Curve {52} = Num_dom_right_up Using Progression A*Prog_dom_right_up;
Line(53) = {23,46};    Transfinite Curve {53} = Num_dom_right_up Using Progression A*Prog_dom_right_up;
Line(54) = {35,19};    Transfinite Curve {54} = Num_dom_right_up Using Progression A*Prog_dom_right_up;

Line(55) = {21,41};    Transfinite Curve {55} = Num_dom_left_up Using Progression B*Prog_dom_left_up;
Line(56) = {40,33};    Transfinite Curve {56} = Num_dom_left_up Using Progression B*Prog_dom_left_up;
Line(57) = {39,32};    Transfinite Curve {57} = Num_dom_left_up Using Progression B*Prog_dom_left_up;
Line(58) = {20,38};    Transfinite Curve {58} = Num_dom_left_up Using Progression B*Prog_dom_left_up;

Line(59) = {20,39};    Transfinite Curve {59} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(60) = {38,32};    Transfinite Curve {60} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(61) = {37,31};    Transfinite Curve {61} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(62) = {36,24};    Transfinite Curve {62} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(63) = {35,23};    Transfinite Curve {63} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(64) = {19,46};    Transfinite Curve {64} = Num_dom_left_left Using Progression B*Prog_dom_left_left;

Line(65) = {21,40};    Transfinite Curve {65} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(66) = {41,33};    Transfinite Curve {66} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(67) = {42,34};    Transfinite Curve {67} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(68) = {43,25};    Transfinite Curve {68} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(69) = {44,26};    Transfinite Curve {69} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(70) = {22,45};    Transfinite Curve {70} = Num_dom_left_left Using Progression B*Prog_dom_left_left;

Line(71) = {39,40};    Transfinite Curve {71} = Num_box Using Progression Prog_box;
Line(72) = {37,38};    Transfinite Curve {72} = Num_box Using Progression Prog_box;
Line(73) = {41,42};    Transfinite Curve {73} = Num_box Using Progression Prog_box;
Line(74) = {43,44};    Transfinite Curve {74} = Num_box Using Progression Prog_box;
Line(75) = {45,46};    Transfinite Curve {75} = Num_box Using Progression Prog_box;
Line(76) = {35,36};    Transfinite Curve {76} = Num_box Using Progression Prog_box;


Line(77) = {37,47};    Transfinite Curve {77} = Num_dom_center Using Progression A*Prog_dom_center;
Line(78) = {47,36};    Transfinite Curve {78} = Num_dom_center Using Progression B*Prog_dom_center;
Line(79) = {31,48};    Transfinite Curve {79} = Num_dom_center Using Progression A*Prog_dom_center;
Line(80) = {34,49};    Transfinite Curve {80} = Num_dom_center Using Progression A*Prog_dom_center;
Line(81) = {42,50};    Transfinite Curve {81} = Num_dom_center Using Progression A*Prog_dom_center;
Line(82) = {48,24};    Transfinite Curve {82} = Num_dom_center Using Progression B*Prog_dom_center;
Line(83) = {49,25};    Transfinite Curve {83} = Num_dom_center Using Progression B*Prog_dom_center;
Line(84) = {50,43};    Transfinite Curve {84} = Num_dom_center Using Progression B*Prog_dom_center;

Line(85) = {47,48};    Transfinite Curve {85} = Num_dom_left_left Using Progression B*Prog_dom_left_left;
Line(86) = {49,50};    Transfinite Curve {86} = Num_dom_left_left Using Progression A*Prog_dom_left_left;

Line(87) = {48,49};    Transfinite Curve {87} = Num_box Using Progression Prog_box;


// ---- SURFACES AND LINES LOOPS ----   made on the GUI


Curve Loop(1) = {-54, -64, 53, 63};

Plane Surface(1) = {1};
Curve Loop(2) = {-28, -53, 75, 52};

Plane Surface(2) = {2};
Curve Loop(3) = {-69, 51, 70, -52};

Plane Surface(3) = {3};
Curve Loop(4) = {76, 62, -25, -63};

Plane Surface(4) = {4};
Curve Loop(5) = {-27, 69, 74, -68};

Plane Surface(5) = {5};
Curve Loop(6) = {-4, -9, 8, 12};

Plane Surface(6) = {6};
//+
Curve Loop(7) = {-3, -12, 7, 11};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {6, -11, -2, 10};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {-10, 5, 9, -1};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {-30, 25, 29, -5};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {28, -29, -8, 32};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {26, -31, -6, 30};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {31, 27, -32, -7};
//+
Plane Surface(13) = {13};
//+
Curve Loop(17) = {60, -41, -61, 72};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {-43, 67, 73, -66};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {-14, -23, 18, 22};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {23, 19, -24, -15};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {-16, -21, 20, 24};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {17, -22, -13, 21};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {-38, 33, 37, -17};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {34, -39, -18, 38};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {35, -40, -19, 39};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {-20, -37, 36, 40};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {41, -46, -33, 45};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {46, 42, -47, -34};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {47, 43, -48, -35};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {-45, 44, 48, -36};
//+
Plane Surface(30) = {30};
//+
Curve Loop(31) = {-58, -60, 57, 59};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {71, 56, -42, -57};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {-65, 55, 66, -56};
//+
Plane Surface(33) = {33};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {11};
//+
Transfinite Surface {10};
//+
Transfinite Surface {12};
//+
Transfinite Surface {13};
//+
Transfinite Surface {7};
//+
Transfinite Surface {6};
//+
Transfinite Surface {9};
//+
Transfinite Surface {8};
//+

//+
Transfinite Surface {17};
//+
Transfinite Surface {18};
//+
Transfinite Surface {29};
//+
Transfinite Surface {30};
//+
Transfinite Surface {27};
//+
Transfinite Surface {28};
//+
Transfinite Surface {25};
//+
Transfinite Surface {26};
//+
Transfinite Surface {23};
//+
Transfinite Surface {24};
//+
Transfinite Surface {20};
//+
Transfinite Surface {21};
//+
Transfinite Surface {22};
//+
Transfinite Surface {19};
//+
Transfinite Surface {31};
//+
Transfinite Surface {32};
//+
Transfinite Surface {33};
//+
Curve Loop(34) = {-62, 82, 85, -78};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {-77, -85, 79, 61};
//+
Plane Surface(35) = {35};
//+
Curve Loop(36) = {-44, -79, -87, 80};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {-82, -26, 83, 87};
//+
Plane Surface(37) = {37};
//+
Curve Loop(38) = {-67, -80, -86, 81};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {-83, 68, 84, 86};
//+
Plane Surface(39) = {39};
//+
Transfinite Surface {34};
//+
Transfinite Surface {35};
//+
Transfinite Surface {36};
//+
Transfinite Surface {37};
Transfinite Surface {38};
Transfinite Surface {39};
Recombine Surface "*";

// PHYSICAL PARAMETERS

Physical Curve(1) = {59,71,65};   //inlet
Physical Curve(2) = {64,75,70};  //outlet 
Physical Curve(3) = {13,14,15,16}; //emitter
Physical Curve(4) = {1,2,3,4};//collector

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Surface(7) = {7};
Physical Surface(8) = {8};
Physical Surface(9) = {9};
Physical Surface(10) = {10};
Physical Surface(11) = {11};
Physical Surface(12) = {12};
Physical Surface(13) = {13};

Physical Surface(17) = {17};
Physical Surface(18) = {18};
Physical Surface(19) = {19};
Physical Surface(20) = {20};
Physical Surface(21) = {21};
Physical Surface(22) = {22};
Physical Surface(23) = {23};
Physical Surface(24) = {24};
Physical Surface(25) = {25};
Physical Surface(26) = {26};
Physical Surface(27) = {27};
Physical Surface(28) = {28};
Physical Surface(29) = {29};
Physical Surface(30) = {30};
Physical Surface(31) = {31};
Physical Surface(32) = {32};
Physical Surface(33) = {33};
Physical Surface(34) = {34};
Physical Surface(35) = {35};
Physical Surface(36) = {36};
Physical Surface(37) = {37};
Physical Surface(38) = {38};
Physical Surface(39) = {39};

