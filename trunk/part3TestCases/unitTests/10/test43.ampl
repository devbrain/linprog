var x1 >=0 ;
var x2 >=0 ;
var x3 >=0 ;
var x4 >=0 ;
var x5 >=0 ;
var x6 >=0 ;
var x7 >=0 ;
var x8 >=0 ;
var x9 >=0 ;
var x10 >=0 ;
var x11 >=0 ;
var x12 >=0 ;
var x13 >=0 ;
maximize obj: 0.0  -3.0 * x1   + 3.0 * x2   -2.0 * x3   + 1.0 * x4   -2.0 * x5   -2.0 * x6 ;
c1: x7 = 21.0  -3.0 * x1  -6.0 * x2  -7.0 * x3  -3.0 * x4  + 9.0 * x5  -1.0 * x6 ;
c2: x8 = 40.0  + 1.0 * x1  -8.0 * x2  + 1.0 * x3  -10.0 * x4  -2.0 * x5  + 8.0 * x6 ;
c3: x9 = 71.0  -5.0 * x1  -7.0 * x2  -9.0 * x3  -7.0 * x4  + 1.0 * x5  -9.0 * x6 ;
c4: x10 = -4.0  + 7.0 * x1  -8.0 * x2  + 1.0 * x3  + 3.0 * x4  + 10.0 * x5  + 7.0 * x6 ;
c5: x11 = -2.0  -4.0 * x1  + 2.0 * x2  -10.0 * x3  -2.0 * x4  + 10.0 * x5  + 10.0 * x6 ;
c6: x12 = 66.0  -8.0 * x1  -9.0 * x2  -5.0 * x3  + 0.0 * x4  + 0.0 * x5  -10.0 * x6 ;
c7: x13 = -21.0  + 7.0 * x1  -9.0 * x2  + 10.0 * x3  + 4.0 * x4  + 9.0 * x5  + 8.0 * x6 ;
solve; 
display 0.0  -3.0 * x1   + 3.0 * x2   -2.0 * x3   + 1.0 * x4   -2.0 * x5   -2.0 * x6 ;
 
 end; 