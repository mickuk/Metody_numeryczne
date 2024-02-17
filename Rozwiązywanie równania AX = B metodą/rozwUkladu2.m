function X = rozwUkladu2(A, B)


% funkcja rozwiązuje układ równań metodą Choleskiego-Banasiewicza poprzez
% stosowanie funkcji rref

[L, D, LT, info] = rozkladChol(A);


% LZ = B  => wyznaczamy Z (czyli D(L^T)X) --> inaczej Z = L^(-1)*B
% następnie
% DY = Z  => wyznaczamy Y 
% i na końcu
% (L^T)X = Y  => wyznaczamy X


A1 = [L B];
Z = rref(A1);
Z = Z(:, size(A1,2) - size(B,2)+1:size(A1,2));

A2 = [D Z];
Y = rref(A2);
Y = Y(:, size(A2,2) - size(Z,2)+1: size(A2,2));

A3 = [LT Y];
X = rref(A3);
X = X(:, size(A3,2)- size(Y,2)+1: size(A3,2));
