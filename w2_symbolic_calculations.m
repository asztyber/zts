syms s;
d = 0;

A = [0 1; -1 -d];
B = [0;1];
C = [1 0];
D = [0];

G = C*inv(eye(2)*s - A)*B + D;

solve(1/G == 0) % zera mianownika transmitancji

eigs(A) % wartoœci w³asne macierzy AG