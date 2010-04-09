//
mode(-1);
lines(0);

c=[90   0   0 ; 1   9   0];    
disp(c, "c =");

disp("[NI,A,Rej,P,R]=confmatrix2ni_mi(c)");
[NI,A,Rej,P,R]=confmatrix2ni_mi(c);
disp(NI, "NI =", A, "A = ", Rej, "Rej = ", P, "P = ", R, "R = ")

