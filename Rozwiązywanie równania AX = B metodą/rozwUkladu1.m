function X = rozwUkladu1(A, B)

% funkcja rozwiązuje układ równań AX = B, przez rozkład Choleskiego
% Banachiewicza poprzez obliczanie odwrotności macierzy



[L, D, LT, info] = rozkladChol(A);

% U nas A = LDL^T, zatem mamy równanie
% LD(L^T)X = B
% (L^T)X = Y
% LDY = B
% DY = Z
% LZ = B

% Będziemy postępować w odwrotnej kolejności, zatem najpierw rozwiązujemy:
% LZ = B  => wyznaczamy Z (czyli D(L^T)X) --> inaczej Z = L^(-1)*B
% następnie
% DY = Z  => wyznaczamy Y 
% i na końcu
% (L^T)X = Y  => wyznaczamy X

Z = L\B;
Y = D\Z;
X = LT\Y;
