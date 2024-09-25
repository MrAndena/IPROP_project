// ===========================================
// ==================================MESH FILE
// ===========================================

// the origin of the mesh is the Point 1


Include "../build/Data4Geo.txt";
Include "../build/NACA_points.txt";

//PARAMETERS

// you are working with a NACA0012 airfoil

r_emi= emitter_radius; //[m]   radius of the emitter 

dist_emi_col = distance_emitter_collector;            //  [m] distance emitter collector (starting from the perimiter)
dist_emi_up_down = distance_emitter_up_bottom;        //  [m] distance between the surface of the emitter and the upper/bottom part of the domain 
dist_Tedge_outlet = distance_trailing_edge_outlet;    //  [m] distance between the trailing edge of the collector and the outlet 
dist_emi_inlet = distance_emitter_inlet;             //   [m] distance between the surface of the emitter and the inlet


H =350*r_emi;           // half of the height of the box that wraps the airfoil
W = 550*r_emi;          // half of the height at which accour a separation of the domain
base = 0.03*chord_length;           // the quantity to add at naca chord length to find the length of the rectangle

L_emi_in = 50*r_emi;       // [m] half the length of the inner square that wraps the emitter
L_emi_out = 150*r_emi;      // [m] half the length of the outer square that wraps the emitter
                          // NB: for what concern the emitter, there is a second other box that wraps it; it has the dimension of L_col_out

mesh_ref_1 = 1;           // dummy mesh ref

// transfinite parameters
A =1.05;
B = 0.95;

N_A = 51;
Prog_A = 1;

N_B = 46;
Prog_B = 1;

N_C = 26;
Prog_C = 1;

N_D = 11;
Prog_D = 1;

N_E = 21;
Prog_E = 1;

N_F = 26;
Prog_F = 1;

N_G = 16;
Prog_G = 1;

N_I = 16;
Prog_I = 1;

N_L = 31;
Prog_L = 1;

N_M = 26;
Prog_M = 1;

N_N = 16;
Prog_N =1;

N_O = 26;
Prog_O =1;




Point(256) = {-base, H, 0, mesh_ref_1};
Point(257) = {-base, -H, 0, mesh_ref_1};
Point(258) = { chord_length, -H, 0, mesh_ref_1};
Point(259) = { chord_length, H, 0, mesh_ref_1};

//EMITTER POINTS
Point(239) = {-dist_emi_col - r_emi, 0, 0, mesh_ref_1}; // center emitter
Point(240) = {-dist_emi_col - r_emi + r_emi/Sqrt(2), r_emi/Sqrt(2), 0, mesh_ref_1};
Point(241) = {-dist_emi_col - r_emi - r_emi/Sqrt(2), r_emi/Sqrt(2), 0, mesh_ref_1};
Point(242) = {-dist_emi_col - r_emi - r_emi/Sqrt(2), -r_emi/Sqrt(2), 0, mesh_ref_1};
Point(243) = {-dist_emi_col - r_emi + r_emi/Sqrt(2), -r_emi/Sqrt(2), 0, mesh_ref_1};

Point(244) = {-dist_emi_col - r_emi + L_emi_in, L_emi_in, 0, mesh_ref_1};
Point(245) = {-dist_emi_col - r_emi - L_emi_in, L_emi_in, 0, mesh_ref_1};
Point(246) = {-dist_emi_col - r_emi - L_emi_in, -L_emi_in, 0, mesh_ref_1};
Point(247) = {-dist_emi_col - r_emi + L_emi_in, -L_emi_in, 0, mesh_ref_1};

Point(248) = {-dist_emi_col - r_emi + L_emi_out, L_emi_out, 0, mesh_ref_1};
Point(249) = {-dist_emi_col - r_emi - L_emi_out, L_emi_out, 0, mesh_ref_1};
Point(250) = {-dist_emi_col - r_emi - L_emi_out, -L_emi_out, 0, mesh_ref_1};
Point(251) = {-dist_emi_col - r_emi + L_emi_out, -L_emi_out, 0, mesh_ref_1};

Point(252) = {-dist_emi_col - r_emi + H, H, 0, mesh_ref_1};
Point(253) = {-dist_emi_col - r_emi - H, H, 0, mesh_ref_1};
Point(254) = {-dist_emi_col - r_emi - H, -H, 0, mesh_ref_1}; 
Point(255) = {-dist_emi_col - r_emi + H, -H, 0, mesh_ref_1};


//RECTANGULAR DOMAIN POINTS
Point(294) = {chord_length + dist_Tedge_outlet, 0, 0, mesh_ref_1};

Point(260) = {chord_length + dist_Tedge_outlet,r_emi+dist_emi_up_down, 0, mesh_ref_1};
Point(261) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, r_emi+dist_emi_up_down, 0, mesh_ref_1};
Point(262) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, -r_emi-dist_emi_up_down, 0, mesh_ref_1};
Point(263) = {chord_length + dist_Tedge_outlet,-r_emi-dist_emi_up_down, 0, mesh_ref_1};

Point(264) = {chord_length + dist_Tedge_outlet, H, 0, mesh_ref_1};
Point(265) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, H, 0, mesh_ref_1};
Point(266) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, -H, 0, mesh_ref_1};
Point(267) = {chord_length + dist_Tedge_outlet,-H, 0, mesh_ref_1};

Point(268) = {chord_length + dist_Tedge_outlet, W, 0, mesh_ref_1};
Point(269) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, W, 0, mesh_ref_1};
Point(270) = {-dist_emi_col - 2*r_emi - dist_emi_inlet, -W, 0, mesh_ref_1};
Point(271) = {chord_length + dist_Tedge_outlet,-W, 0, mesh_ref_1};

Point(272) = {chord_length, W, 0, mesh_ref_1};
Point(273) = {chord_length, -W, 0, mesh_ref_1};
Point(274) = {chord_length, r_emi + dist_emi_up_down, 0, mesh_ref_1};
Point(275) = {chord_length, -r_emi - dist_emi_up_down, 0, mesh_ref_1};

Point(276) = {-base, W, 0, mesh_ref_1};
Point(277) = {-base, -W, 0, mesh_ref_1};
Point(278) = {-base, r_emi + dist_emi_up_down, 0, mesh_ref_1};
Point(279) = {-base, -r_emi - dist_emi_up_down, 0, mesh_ref_1};

Point(280) = {-0.4*dist_emi_col, H, 0, mesh_ref_1};
Point(281) = {-0.4*dist_emi_col, -H, 0, mesh_ref_1};
Point(282) = {-0.4*dist_emi_col, W, 0, mesh_ref_1};
Point(283) = {-0.4*dist_emi_col, -W, 0, mesh_ref_1};
Point(284) = {-0.4*dist_emi_col, r_emi + dist_emi_up_down, 0, mesh_ref_1};
Point(285) = {-0.4*dist_emi_col, -r_emi - dist_emi_up_down, 0, mesh_ref_1};

Point(286) = {-dist_emi_col - r_emi + H, W, 0, mesh_ref_1};
Point(287) = {-dist_emi_col - r_emi + H, -W, 0, mesh_ref_1};
Point(288) = {-dist_emi_col - r_emi + H, r_emi + dist_emi_up_down, 0, mesh_ref_1};
Point(289) = {-dist_emi_col - r_emi + H, -r_emi - dist_emi_up_down, 0, mesh_ref_1};

Point(290) = {-dist_emi_col - r_emi - H, r_emi + dist_emi_up_down, 0, mesh_ref_1};
Point(291) = {-dist_emi_col - r_emi - H, -r_emi - dist_emi_up_down, 0, mesh_ref_1};
Point(292) = {-dist_emi_col - r_emi - H, W, 0, mesh_ref_1};
Point(293) = {-dist_emi_col - r_emi - H, -W, 0, mesh_ref_1};

//DOMAIN LINES
 
// vertical lines from left to right
Line(9) = {261, 269};
Line(10) = {269, 265};
Line(11) = {265, 266};
Line(12) = {266, 270};
Line(13) = {270, 262};
Line(14) = {290, 292};
Line(15) = {292, 253};
Line(17) = {253, 254};
Line(18) = {254, 293};
Line(19) = {293, 291};
Line(20) = {288, 286};
Line(21) = {286, 252};
Line(22) = {252, 255};
Line(23) = {255, 287};
Line(24) = {287, 289};
Line(25) = {284, 282};
Line(26) = {282, 280};
Line(27) = {280, 281};
Line(28) = {281, 283};
Line(29) = {283, 285};
Line(30) = {278, 276};
Line(31) = {276, 256};
Line(33) = {256, 257};
Line(34) = {257, 277};
Line(35) = {277, 279};
Line(36) = {274, 272};
Line(37) = {272, 259};
Line(38) = {259, 120};
Line(39) = {120, 258};
Line(40) = {258, 273};
Line(41) = {273, 275};
Line(42) = {260, 268};
Line(43) = {268, 264};
Line(45) = {267, 271};
Line(46) = {271, 263};
Line(83) = {120, 294};
Line(84) = {264, 294};
Line(85) = {294, 267};

// orizzontal lines top to bottom

Line(47) = {261, 290};
Line(48) = {290, 288};
Line(49) = {288, 284};
Line(50) = {284, 278};
Line(51) = {278, 274};
Line(52) = {274, 260};
Line(53) = {269, 292};
Line(54) = {292, 286};
Line(55) = {286, 282};
Line(56) = {282, 276};
Line(57) = {276, 272};
Line(58) = {272, 268};
Line(59) = {265, 253};
Line(60) = {253, 252};
Line(61) = {252, 280};
Line(62) = {280, 256};
Line(63) = {256, 259};
Line(64) = {259, 264};
Line(65) = {266, 254};
Line(66) = {254, 255};
Line(67) = {255, 281};
Line(68) = {281, 257};
Line(69) = {257, 258};
Line(70) = {258, 267};
Line(71) = {270, 293};
Line(72) = {293, 287};
Line(73) = {287, 283};
Line(74) = {283, 277};
Line(75) = {277, 273};
Line(76) = {273, 271};
Line(77) = {262, 291};
Line(78) = {291, 289};
Line(79) = {289, 285};
Line(80) = {285, 279};
Line(81) = {279, 275};
Line(82) = {275, 263};

//CIRCULAR ARCS
Circle(5) = {240   , 239   , 241};
Circle(6) = {241   , 239   , 242};
Circle(7) = {242   , 239   , 243};
Circle(8) = {243   , 239   , 240};

Line(86) = {247, 244};
Line(87) = {244, 245};
Line(88) = {245, 246};
Line(89) = {246, 247};
Line(90) = {240, 244};
Line(91) = {241, 245};
Line(92) = {242, 246};
Line(93) = {243, 247};
Line(94) = {248, 249};
Line(95) = {249, 250};
Line(96) = {250, 251};
Line(97) = {251, 248};
Line(98) = {244, 248};
Line(99) = {245, 249};
Line(100) = {246, 250};
Line(101) = {247, 251};
Line(102) = {248, 252};
Line(103) = {249, 253};
Line(104) = {250, 254};
Line(105) = {251, 255};

//AIRFOIL CURVE
Spline(1) = {1:20};
Spline(2) = {20:120};
Spline(3) = {120:220};
Spline(4) = {220:238,1};

Line(106) = {20,256};
Line(107) = {220, 257};

// TRANSFINITE PROCESS

// most up and bottom domain vertical lines
Transfinite Curve {9} =  N_D  Using Progression Prog_D;
Transfinite Curve {14} = N_D  Using Progression Prog_D;
Transfinite Curve {20} = N_D  Using Progression Prog_D;
Transfinite Curve {25} = N_D  Using Progression Prog_D;
Transfinite Curve {30} = N_D  Using Progression Prog_D;
Transfinite Curve {36} = N_D  Using Progression Prog_D;
Transfinite Curve {42} = N_D  Using Progression Prog_D;
Transfinite Curve {13} = N_D  Using Progression Prog_D;
Transfinite Curve {19} = N_D  Using Progression Prog_D;
Transfinite Curve {24} = N_D  Using Progression Prog_D;
Transfinite Curve {29} = N_D  Using Progression Prog_D;
Transfinite Curve {35} = N_D  Using Progression Prog_D;
Transfinite Curve {41} = N_D  Using Progression Prog_D;
Transfinite Curve {46} = N_D  Using Progression Prog_D;


// middle up bottom domain vertical lines
Transfinite Curve {10} = N_E Using Progression Prog_E;
Transfinite Curve {15} = N_E Using Progression Prog_E;
Transfinite Curve {21} = N_E Using Progression Prog_E;
Transfinite Curve {26} = N_E Using Progression Prog_E;
Transfinite Curve {31} = N_E Using Progression Prog_E;
Transfinite Curve {37} = N_E Using Progression Prog_E;
Transfinite Curve {43} = N_E Using Progression Prog_E;
Transfinite Curve {12} = N_E Using Progression Prog_E;
Transfinite Curve {18} = N_E Using Progression Prog_E;
Transfinite Curve {23} = N_E Using Progression Prog_E;
Transfinite Curve {28} = N_E Using Progression Prog_E;
Transfinite Curve {34} = N_E Using Progression Prog_E;
Transfinite Curve {40} = N_E Using Progression Prog_E;
Transfinite Curve {45} = N_E Using Progression Prog_E;

// left left domain orizzontal lines
Transfinite Curve {47} = N_F Using Progression Prog_F;
Transfinite Curve {53} = N_F Using Progression Prog_F;
Transfinite Curve {59} = N_F Using Progression Prog_F;
Transfinite Curve {65} = N_F Using Progression Prog_F;
Transfinite Curve {71} = N_F Using Progression Prog_F;
Transfinite Curve {77} = N_F Using Progression Prog_F;


// left middle domain orizzontal lines
Transfinite Curve {48} = N_C Using Progression Prog_C;
Transfinite Curve {54} = N_C Using Progression Prog_C;
Transfinite Curve {60} = N_C Using Progression Prog_C;
Transfinite Curve {66} = N_C Using Progression Prog_C;
Transfinite Curve {72} = N_C Using Progression Prog_C;
Transfinite Curve {78} = N_C Using Progression Prog_C;

// middle domain orizzontal lines
// left part
Transfinite Curve {49} = N_G Using Progression Prog_G;
Transfinite Curve {55} = N_G Using Progression Prog_G;
Transfinite Curve {61} = N_G Using Progression Prog_G;
Transfinite Curve {67} = N_G Using Progression Prog_G;
Transfinite Curve {73} = N_G Using Progression Prog_G;
Transfinite Curve {79} = N_G Using Progression Prog_G;

//right part
Transfinite Curve {50} = N_I Using Progression Prog_I;
Transfinite Curve {56} = N_I Using Progression Prog_I;
Transfinite Curve {62} = N_I Using Progression Prog_I;
Transfinite Curve {68} = N_I Using Progression Prog_I;
Transfinite Curve {74} = N_I Using Progression Prog_I;
Transfinite Curve {80} = N_I Using Progression Prog_I;

// orizzontal lines over below airfoil
Transfinite Curve {51} = N_B Using Progression Prog_B;
Transfinite Curve {57} = N_B Using Progression Prog_B;
Transfinite Curve {63} = N_B Using Progression Prog_B;
Transfinite Curve {69} = N_B Using Progression Prog_B;
Transfinite Curve {75} = N_B Using Progression Prog_B;
Transfinite Curve {81} = N_B Using Progression Prog_B;

// airfoil up dpwn
Transfinite Curve {2} = N_B Using Progression Prog_B;
Transfinite Curve {3} = N_B Using Progression Prog_B;

// domain right up and bottom orizzontal lines
Transfinite Curve {52} = N_L Using Progression Prog_L;
Transfinite Curve {58} = N_L Using Progression Prog_L;
Transfinite Curve {64} = N_L Using Progression Prog_L;
Transfinite Curve {83} = N_L Using Progression Prog_L;
Transfinite Curve {70} = N_L Using Progression Prog_L;
Transfinite Curve {76} = N_L Using Progression Prog_L;
Transfinite Curve {82} = N_L Using Progression Prog_L;

// arfoil part obliquos edge & co
Transfinite Curve {84} = N_A Using Progression B*Prog_A;
Transfinite Curve {38} = N_A Using Progression B*Prog_A;
Transfinite Curve {106} = N_A Using Progression A*Prog_A;
Transfinite Curve {107} = N_A Using Progression A*Prog_A;
Transfinite Curve {39} = N_A Using Progression A*Prog_A;
Transfinite Curve {85} = N_A Using Progression A*Prog_A;

// THE LEFT AIRFOIL EDGE
Transfinite Curve {1} = 14 Using Progression Prog_C;
Transfinite Curve {4} = 13 Using Progression Prog_C;
Transfinite Curve {33} = N_C Using Progression Prog_C;

// last edges around emitter
Transfinite Curve {27} = N_C Using Progression Prog_C;
Transfinite Curve {22} = N_C Using Progression Prog_C;
Transfinite Curve {17} = N_C Using Progression Prog_C;
Transfinite Curve {11} = N_C Using Progression Prog_C;

Transfinite Curve {94} = N_C Using Progression Prog_C;
Transfinite Curve {95} = N_C Using Progression Prog_C;
Transfinite Curve {97} = N_C Using Progression Prog_C;
Transfinite Curve {96} = N_C Using Progression Prog_C;

Transfinite Curve {88} = N_C Using Progression Prog_C;
Transfinite Curve {87} = N_C Using Progression Prog_C;
Transfinite Curve {86} = N_C Using Progression Prog_C;
Transfinite Curve {89} = N_C Using Progression Prog_C;

Transfinite Curve {6} = N_C Using Progression Prog_C;
Transfinite Curve {5} = N_C Using Progression Prog_C;
Transfinite Curve {8} = N_C Using Progression Prog_C;
Transfinite Curve {7} = N_C Using Progression Prog_C;

Transfinite Curve {103} = N_M Using Progression Prog_M;
Transfinite Curve {102} = N_M Using Progression Prog_M;
Transfinite Curve {105} = N_M Using Progression Prog_M;
Transfinite Curve {104} = N_M Using Progression Prog_M;

Transfinite Curve {99} = N_N Using Progression Prog_N;
Transfinite Curve {98} = N_N Using Progression Prog_N;
Transfinite Curve {101} = N_N Using Progression Prog_N;
Transfinite Curve {100} = N_N Using Progression Prog_N;

Transfinite Curve {91} = N_O Using Progression A*Prog_O;
Transfinite Curve {90} = N_O Using Progression A*Prog_O;
Transfinite Curve {93} = N_O Using Progression A*Prog_O;
Transfinite Curve {92} = N_O Using Progression A*Prog_O;
//+
Curve Loop(1) = {9, 53, -14, -47};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {14, 54, -20, -48};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {20, 55, -25, -49};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {25, 56, -30, -50};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {30, 57, -36, -51};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {36, 58, -42, -52};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {10, 59, -15, -53};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {15, 60, -21, -54};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {21, 61, -26, -55};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {26, 62, -31, -56};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {31, 63, -37, -57};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {37, 64, -43, -58};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {11, 65, -17, -59};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {17, -104, -95, 103};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {-103, -60, 102, -94};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {-97, -102, -22, 105};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {66, -105, -96, 104};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {-100, 95, 99, -88};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {-99, 94, 98, -87};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {-86, -98, 97, 101};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {100, 96, -101, -89};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {88, -92, -6, 91};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {-91, 87, 90, -5};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {-8, -90, 86, 93};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {92, 89, -93, -7};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {22, 67, -27, -61};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {27, 68, -33, -62};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {33, -107, 4, 1, 106};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {-106, -63, -38, 2};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {69, -39, 3, 107};
//+
Plane Surface(30) = {30};
//+
Curve Loop(31) = {38, 83, -84, -64};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {39, 70, -85, -83};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {12, 71, -18, -65};
//+
Plane Surface(33) = {33};
//+
Curve Loop(34) = {18, 72, -23, -66};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {23, 73, -28, -67};
//+
Plane Surface(35) = {35};
//+
Curve Loop(36) = {28, 74, -34, -68};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {34, 75, -40, -69};
//+
Plane Surface(37) = {37};
//+
Curve Loop(38) = {40, 76, -45, -70};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {13, 77, -19, -71};
//+
Plane Surface(39) = {39};
//+
Curve Loop(40) = {19, 78, -24, -72};
//+
Plane Surface(40) = {40};
//+
Curve Loop(41) = {24, 79, -29, -73};
//+
Plane Surface(41) = {41};
//+
Curve Loop(42) = {29, 80, -35, -74};
//+
Plane Surface(42) = {42};
//+
Curve Loop(43) = {35, 81, -41, -75};
//+
Plane Surface(43) = {43};
//+
Curve Loop(44) = {41, 82, -46, -76};
Plane Surface(44) = {44};
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
Transfinite Surface {6};
//+
Transfinite Surface {7};
//+
Transfinite Surface {8};
//+
Transfinite Surface {9};
//+
Transfinite Surface {10};
//+
Transfinite Surface {11};
//+
Transfinite Surface {12};
//+
Transfinite Surface {13};
//+
Transfinite Surface {14};
//+
Transfinite Surface {15};
//+
Transfinite Surface {16};
//+
Transfinite Surface {17};
//+
Transfinite Surface {18};
//+
Transfinite Surface {19};
//+
Transfinite Surface {20};
//+
Transfinite Surface {21};
//+
Transfinite Surface {22};
//+
Transfinite Surface {23};
//+
Transfinite Surface {24};
//+
Transfinite Surface {25};
//+
Transfinite Surface {26};
//+
Transfinite Surface {27};
//+
Transfinite Surface {28}={257,256,20,220};
//+
Transfinite Surface {29};
//+
Transfinite Surface {30};
//+
Transfinite Surface {31};
//+
Transfinite Surface {32};
//+
Transfinite Surface {33};
//+
Transfinite Surface {34};
//+
Transfinite Surface {35};
//+
Transfinite Surface {36};
//+
Transfinite Surface {37};
//+
Transfinite Surface {38};
//+
Transfinite Surface {39};
//+
Transfinite Surface {40};
//+
Transfinite Surface {41};
//+
Transfinite Surface {42};
//+
Transfinite Surface {43};
//+
Transfinite Surface {44};
Recombine Surface "*";

// PHYSICAL PARAMETERS

Physical Curve(1) = {9,10,11,12,13};   //inlet
Physical Curve(2) = {42,43,84,85,45,46};  //outlet 
Physical Curve(3) = {5,6,7,8}; //emitter
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
Physical Surface(14) = {14};
Physical Surface(15) = {15};
Physical Surface(16) = {16};
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
Physical Surface(40) = {40};
Physical Surface(41) = {41};
Physical Surface(42) = {42};
Physical Surface(43) = {43};
Physical Surface(44) = {44};
