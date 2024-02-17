function [L, D, LT, info] = rozkladChol(A)
% rozkład LDL^(T) macierzy A
% rozkład istnieje, gdy info = 1, wpp info = 0


[m,n] = size(A);
info = 0;
% jeżeli macierz A spełnia warunki zadania, to pod koniec info = 1
% sprawdzamy, czy macierz jest kwadratowa
if m ~= n
    return
end

% Na początek przyjmiemy, że nasze macierze (w tym macierz C - potrzebna
% później przy rozkładzie LDL^T) są macierzami jednostkowymi.
% W macierzy L będziemy w algorytmie zmieniać wyrazy tylko pod dolną
% przekątną (przekątna to jedynki, a nad przekątną zera), a w macierzy D
% (diagonalnej) tylko wyrazy na przekątnej. Macierz LT to transponowana L
L = eye(n);
D = eye(n);
C = eye(n);
LT = eye(n);

% sprawdzamy, czy macierz jest symetryczna
%for i = 1:m
%    for j = 1:i
%        if A(i,j) ~= A(j,i)
%            return
%        end
%    end
%end

% powyższy kod nie wykorzystuje wektoryzacji

for i = 1:m
    if A(i,:) ~= A(:,i)'
        return
    end
end

% sprawdzamy, czy macierz jest dodatnio określona
for i = 1:m
    if det(A(1:i,1:i)) <= 0
        return
    end
end


% i tutaj implementujemy algorytm (wzięty z książki o metodach
% numerycznych)

D(1,1) = A(1,1);
for i = 2:n
    for j = 1:i-1
        suma1 = 0;
        for k = 1:j-1
            if j == 1
                continue
            else
            suma1 = suma1 + C(i,k)*L(j,k);
            end
        end
        L(i,j) = (A(i,j) - suma1)/D(j,j);
        C(i,j) = D(j,j)*L(i,j);
    end
    suma2=0;
    for k = 1:i-1
        suma2 = suma2 + C(i,k)*L(i,k);
    end
    D(i,i) = A(i,i) - suma2;
end


info = 1;  % jeżeli dotarliśmy do tego momentu, to metoda działa
LT = L';

end




