function [x] = solve_chol(A, B)
% Funkcja ta zwraca nam X, które jest rozwiązaniem równania LDL^TX = B,
% gdzie A = LDL^T

[L,D, LT, info] = rozkladChol(A);

n = size(L,1);
m = size(B,2);
x = zeros(n,m);

y = zeros(n,m);
d = zeros(n,1);
for i = 1:n
    d(i,1) = D(i,i);
end


for i = 1:n                         % wyliczamy i-tą współrzędną w i-tym wierszu na podstawie współrzędnych od 1 do i-1
    sum = L(i, 1:i-1)*y(1:i-1,:);   % iloczyn wiersza i-tego aż do przekątnej oraz
    y(i,:) = B(i,:)-sum;            % rozwiązań układu aż do i-tego elementu
end
L = L';
y = y./d;                       % rozwiązanie D*x=y
for i = flip(1:n)               % analogicznie jak wcześniej tylko dla L^T
    sum = L(i, i+1:n)*x(i+1:n,:);
    x(i,:) = y(i,:)-sum;
end
L = L';
y = zeros(n,m);

end