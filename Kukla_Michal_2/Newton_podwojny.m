function [X, info, k] = Newton_podwojny(f, x0, a, b)

% funkcja wyznacza podwójny pierwiastek funkcji najbliżej położenia
% początkowego

% wszystko analogicznie jak w Newton_pojedynczy

eps = 1e-10;
info = 0;
X = 100;
k = 1;

lista = zeros(1,1000);

lista(1) = x0;

if f(a)*f(b) >= 0
    return
end

for k = 1:1000
    if pochodna(f, lista(k)) == 0
        X = lista(k);
        return
    end
    lista(k+1) = lista(k) - 2*f(lista(k))/pochodna(f, lista(k));
    if abs(lista(k+1) - lista(k)) < eps
        X = lista(k+1);
        info = 1;
        return
    end
end

end