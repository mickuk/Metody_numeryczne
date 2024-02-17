function [A,B] = gen_mac_sym_dod(n,m)
% funkcja generuje macierz symetryczną dodatnio określoną i macierz B
% funkcja przyjmuje argument n (ile wierszy chcemy dla macierzy A)
% oraz argument m (ile kolumn ma mieć macierz B)

A = rand(n,n);
A = 0.5*(A+A');   % konstruuję symetryczną macierz
A = A + n*eye(n); % korzystam z własności, że macierz silnie dominująca
                 % diagonalnie o elementach na przekątnej dodatnich jest
                 % dodatnio określona
B = 10*rand(n, m); 
end