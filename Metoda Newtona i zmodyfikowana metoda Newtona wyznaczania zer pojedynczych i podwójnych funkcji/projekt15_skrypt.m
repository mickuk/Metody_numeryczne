% skrypt testujący

% definiuję funkcje

format long;

a = @(x) x+4;
b = @(x) x^2 - 7;
c = @(x) x^2 - 4;
d = @(x) x^4 - 4*x^3 + 5*x^2 - 4*x + 4;
e = @(x) x^2*(sin(x) - cos(2*x)+1);
f = @(x) x^3*(x+sin(x^2 -1)-1) -1;
g = @(x) x - (1/3)*x^(1/3)-2;
h = @(x) x^2;
i = @(x) x^4 - 4*x^3 + 5*x^2 - 4*x + 4;
j = @(x) (x-3)^2*(x-1)*x;
k1 = @(x) (sin(x)-2)^3*(x+cos(x))^3*(x-2)^2;
l = @(x) (x-2)^2*(x-1);
m = @(x) (x-3)^2*(x-10);

% sprawdzam, jak działają funkcje
% niektóre funkcje przyjmują wartości tylko dodatnie, stąd metoda Newtona
% nie zadziała, bo nie będzie spełniony warunek f(a)*f(b) < 0

% METODA NEWTONA ZER JEDNOKROTNYCH

[A, info, k] = Newton_pojedynczy(a, -3.5,-5, -3)

[B1, info, k] = Newton_pojedynczy(b, 2, 1,3)

[B2, info, k] = Newton_pojedynczy(b, -2, -3, -1)

[C1, info, k] = Newton_pojedynczy(c, 1, 0, 3)

[C2, info, k] = Newton_pojedynczy(c, -1, -3, 0)

[D, info, k] = Newton_pojedynczy(d, 1.5, 0, 2)

[E, info, k] = Newton_pojedynczy(e, 3, 2, 3.5)

[F1, info, k] = Newton_pojedynczy(f, -0.5, -1, 0)

[F2, info, k] = Newton_pojedynczy(f, 1.5, 1, 2)

[G, info, k] = Newton_pojedynczy(g, 3, 0, 4)

% METODA NEWTONA ZER DWUKROTNYCH

[H, info, k] = Newton_podwojny(h, 1, 0,2)

[I, info, k] = Newton_podwojny(i, 1.5, 0, 2)

[J, info, k] = Newton_podwojny(j, 5, 0.5, 6)

[K, info, k] = Newton_podwojny(k1, 1, -0.9,3)

[L, info, k] = Newton_podwojny(l, 2.5, 0,3)

[M, info, k] = Newton_podwojny(m, 5, 0,11)

% Porównanie czasu dla pojedynczego i podwójnego
disp("gergergre\n\n\n\n/ngreg")

tic; [J, info, k] = Newton_pojedynczy(j, 5, 0.5, 6); toc

tic; [K, info, k] = Newton_pojedynczy(k1, 1, -0.9,3); toc

tic; [L, info, k] = Newton_pojedynczy(l, 2.5, 0,3); toc

tic; [M, info, k] = Newton_pojedynczy(m, 5, 0,11); toc


tic; [J, info, k] = Newton_podwojny(j, 5, 0.5, 6); toc

tic; [K, info, k] = Newton_podwojny(k1, 1, -0.9,3); toc

tic; [L, info, k] = Newton_podwojny(l, 2.5, 0,3); toc

tic; [M, info, k] = Newton_podwojny(m, 5, 0,11); toc






