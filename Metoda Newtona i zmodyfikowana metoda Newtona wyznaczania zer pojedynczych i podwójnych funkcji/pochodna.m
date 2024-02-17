function [dx, info] = pochodna(f, x)


% ustawiam wartość h na bardzo małą

info = 0;
h = 1e-6;

dx = (f(x+h) - f(x))/h;

end


