function [X, info, k] = Newton_pojedynczy(f, x0, a, b)

% funkcja wyznacza pojedynczy pierwiastek funkcji najbliżej położenia
% początkowego

% w metodzie Newtona przyjmuje się następujące założenia funkcji f
% 1 - w [a,b] jest dokładnie jeden pierwiastek - sprawdzam w testowym
% 2 - f(a)*f(b) < 0
% 3 - f' i f'' nie zmieniają znaku w przedziale
% jako że punkt 3 jest trudny do pokazania (wymaga znalenienia pierwiastków
% f' i f'' i naszkicowania wykresu), to w pliku testowym będę podawał
% funkcje ze spełnionym warunkiem 3.

k = 1;
eps = 1e-10;
info = 0;
X = 100;

lista = zeros(1,1000);  % tworzę zbiór kolejnych wyrazów ciągu przybliżeń
lista(1) = x0;

if f(a)*f(b) >= 0
    return
end

% stosuję algorytm do metody Newtona

for k = 1:1000
    if pochodna(f, lista(k)) == 0  % za każdym razem sprawdzam, czy nie
        X = lista(k);              % będę dzielić przez zero
        return
    end
    lista(k+1) = lista(k) - f(lista(k))/pochodna(f, lista(k));
    if abs(lista(k+1) - lista(k)) < eps % jeżeli kolejne wyrazy ciągu są
        X = lista(k+1);                 % blisko siebie, to kończymy
        info = 1;
        return
    end
end

end