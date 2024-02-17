% Skrypt testujący funkcję rozkładu Choleskiego-Banachiewicza LDL^T 
% i rozwiązujący układ równań AX = B (macierze A i B są generowane)

clc
clear

disp("Przykłady wymyślone (niewygenerowane): ")
disp("Pierwszy przykład: ")
A1 = [4,-2,-2; -2,5,-1; -2,-1,6];
B1 = eye(3);
[L1, D1, L1T, info] = rozkladChol(A1)

rozwUkladu1(A1, B1)
rozwUkladu2(A1,B1)
solve_chol(A1,B1)
inv(A1)*B1

tic; rozwUkladu1(A1, B1); toc
tic; rozwUkladu2(A1,B1); toc
tic; solve_chol(A1,B1); toc
tic; inv(A1)*B1; toc


disp("Drugi przykład: ")

A2 = [9,-6, 3; -6,5,1; 3,1,11]
B2 = [-9,1,2,3,4; 5,1,2,3,4; 6,1,2,3,4]
[L2, D2, L2T, info] = rozkladChol(A2)

rozwUkladu1(A2,B2)
rozwUkladu2(A2,B2)
solve_chol(A2,B2)
inv(A2)*B2

tic; rozwUkladu1(A2, B2); toc
tic; rozwUkladu2(A2,B2); toc
tic; solve_chol(A2,B2); toc
tic; inv(A2)*B2; toc

disp("Trzeci przykład")


A3 = [1,-2,3,1; -2,5,-8,1; 3,-8,17,-7; 1,1,-7,18];
B3 = [-1,1; -1,2; 3,3; -4,5];
[L3, D3, L3T, info] = rozkladChol(A3)

rozwUkladu1(A3, B3)
rozwUkladu2(A3,B3)
solve_chol(A3,B3)
inv(A3)*B3


tic; rozwUkladu1(A3, B3); toc
tic; rozwUkladu2(A3,B3); toc
tic; solve_chol(A3,B3); toc
tic; inv(A3)*B3; toc

disp("Czwarty przykład")


A4 = [4,2,-2,4; 2,10,-7,4; -2,-7,6,-5; 4,4,-5,18]
B4 = [-12, 1, 0; 9, 34, 4; -9,0,0; 39,-15,-2]
[L4, D4, L4T, info] = rozkladChol(A4)


rozwUkladu1(A4, B4)
rozwUkladu2(A4,B4)
solve_chol(A4,B4)
inv(A4)*B4

tic; rozwUkladu1(A4, B4); toc
tic; rozwUkladu2(A4,B4); toc
tic; solve_chol(A4,B4); toc
tic; inv(A4)*B4; toc







disp("Przykłady z generowanymi macierzami:")
disp("Piąty przykład: ")

[A5, B5] = gen_mac_sym_dod(10, 3);
B5 = eye(10);  % użyjemy tu macierzy jednostkowej

A5
B5

[L5,D5,L5T, info] = rozkladChol(A5)

rozwUkladu1(A5,B5)
rozwUkladu2(A5,B5)
solve_chol(A5,B5)
inv(A5)*B5

tic; rozwUkladu1(A5,B5); toc
tic; rozwUkladu2(A5,B5); toc
tic; solve_chol(A5,B5); toc
tic; inv(A5)*B5; toc

disp("Szósty przykład: ")

[A6, B6] = gen_mac_sym_dod(25,4);

A6
B6

[L6,D6,L6T, info] = rozkladChol(A6)


rozwUkladu1(A6,B6)
rozwUkladu2(A6,B6)
solve_chol(A6,B6)
inv(A6)*B6

tic; rozwUkladu1(A6,B6); toc
tic; rozwUkladu2(A6,B6); toc
tic; solve_chol(A6,B6); toc
tic; inv(A6)*B6; toc

disp("Siódmy przykład: ")

[A7, B7] = gen_mac_sym_dod(50,10);
A7
B7

[L7,D7,L7T, info] = rozkladChol(A7)


rozwUkladu1(A7,B7)
rozwUkladu2(A7,B7)
solve_chol(A7,B7)
inv(A7)*B7

tic; rozwUkladu1(A7,B7); toc
tic; rozwUkladu2(A7,B7); toc
tic; solve_chol(A7,B7); toc
tic; inv(A7)*B7; toc

disp("Ósmy przykład: ")

[A8, B8] = gen_mac_sym_dod(100,100);
A8
B8

[L8,D8,L8T, info] = rozkladChol(A8)


rozwUkladu1(A8,B8)
rozwUkladu2(A8,B8)
solve_chol(A8,B8)
inv(A8)*B8

tic; rozwUkladu1(A8,B8); toc
tic; rozwUkladu2(A8,B8); toc
tic; solve_chol(A8,B8); toc
tic; inv(A8)*B8; toc

% Uznajemy, że wbudowana funkcja inv zwraca poprawną wartość
% Jak widać, przy każdym przykładzie wynik X jest taki sam niezależnie od
% tego, której z czterech metod użyliśmy. Zatem metoda rozwiązywania układu 
% równań przez rozkład LDL^T jest poprawna


% Jednak rozwiązywanie układu równań poprzez najpierw rozkład LDL^T jest
% dłuższe 


% Teraz (zakładając, że wynik funkcji inv jest dokładny) zobaczymy, jak
% bardzo różnią się wyniki funkcji solve_chol() od inv(). Zbadamy tylko 
% solve_chol() (a nie rozwUkladu1() i rozwUkladu2()), ponieważ jest 
% to funkcja napisana samodzielnie, która nie wykorzystuje funkcji "\" oraz
% funkcji rref

disp(" ")
disp("Badanie błędu wynikowego")
A =solve_chol(A1,B1) - inv(A1)*B1
B =solve_chol(A2,B2) - inv(A2)*B2
C =solve_chol(A3,B3) - inv(A3)*B3
D =solve_chol(A4,B4) - inv(A4)*B4
E =solve_chol(A5,B5) - inv(A5)*B5
F =solve_chol(A6,B6) - inv(A6)*B6
G =solve_chol(A7,B7) - inv(A7)*B7
H =solve_chol(A8,B8) - inv(A8)*B8


% W wyniku otrzymaliśmy, że błędy wartości z wyznaczonej macierzy X
% są mniejsze niż 10^(-13) w każdej metodzie.
% Oznacza to, że metoda jest dokładna


sum(A, "all")
sum(B, "all")
sum(C, "all")
sum(D, "all")
sum(E, "all")
sum(F, "all")
sum(G, "all")
sum(H, "all")

% i tutaj badam sumę różnic wszystkich elementów - w najgorszym przypadku
% jest to 10^(-12)



%--------------------------------------------------------------------------
% Teraz ogólne błędy metod

%Wskaźnik uwarunkowania macierzy: 
cond1 =cond(A1)
cond2 =cond(A2)
cond3 =cond(A3)
cond4 =cond(A4)
cond5 =cond(A5)
cond6 =cond(A6)
cond7 =cond(A7)
cond8 =cond(A8)

% jak widać, mimo różnych wielkości macierzy, te prawie najmniejsze były 
% najgorzej uwarunkowane

% Błędy rozkładu

edec1 = norm(A1 -L1*D1*L1T)/norm(A1)
edec2 = norm(A2 -L2*D1*L2T)/norm(A2)
edec3 = norm(A3 -L3*D3*L3T)/norm(A3)
edec4 = norm(A4 -L4*D4*L4T)/norm(A4)
edec5 = norm(A5 -L5*D5*L5T)/norm(A5)
edec6 = norm(A6 -L6*D6*L6T)/norm(A6)
edec7 = norm(A7 -L7*D7*L7T)/norm(A7)
edec8 = norm(A8 -L8*D8*L8T)/norm(A8)

% Błędy generalnie są bardzo małe i wynoszą max 10^-(17)
% Jedynie w przypadku drugiej macierzy było to ponad 2


% Błąd względny

erel1 = norm(solve_chol(A1,B1) - inv(A1)*B1)/norm(inv(A1)*B1)
erel2 = norm(solve_chol(A2,B2) - inv(A2)*B2)/norm(inv(A2)*B2)
erel3 = norm(solve_chol(A3,B3) - inv(A3)*B3)/norm(inv(A3)*B3)
erel4 = norm(solve_chol(A4,B4) - inv(A4)*B4)/norm(inv(A4)*B4)
erel5 = norm(solve_chol(A5,B5) - inv(A5)*B5)/norm(inv(A5)*B5)
erel6 = norm(solve_chol(A6,B6) - inv(A6)*B6)/norm(inv(A6)*B6)
erel7 = norm(solve_chol(A7,B7) - inv(A7)*B7)/norm(inv(A7)*B7)
erel8 = norm(solve_chol(A8,B8) - inv(A8)*B8)/norm(inv(A8)*B8)


% Błąd względny w każdym przypadku jest mały i wynosi max 10^(-14)

% Współczynnik stabilności

wsps1 = norm(solve_chol(A1,B1) - inv(A1)*B1)/(norm(inv(A1)*B1)*cond(A1))
wsps2 = norm(solve_chol(A2,B2) - inv(A2)*B2)/(norm(inv(A2)*B2)*cond(A2))
wsps3 = norm(solve_chol(A3,B3) - inv(A3)*B3)/(norm(inv(A3)*B3)*cond(A3))
wsps4 = norm(solve_chol(A4,B4) - inv(A4)*B4)/(norm(inv(A4)*B4)*cond(A4))
wsps5 = norm(solve_chol(A5,B5) - inv(A5)*B5)/(norm(inv(A5)*B5)*cond(A5))
wsps6 = norm(solve_chol(A6,B6) - inv(A6)*B6)/(norm(inv(A6)*B6)*cond(A6))
wsps7 = norm(solve_chol(A7,B7) - inv(A7)*B7)/(norm(inv(A7)*B7)*cond(A7))
wsps8 = norm(solve_chol(A8,B8) - inv(A8)*B8)/(norm(inv(A8)*B8)*cond(A8))


% maksymalny współczynnik wynosi 10^(-15)

% Współczynnik poprawności

wspp1 = norm(B1 - A1*solve_chol(A1,B1))/(norm(A1)*norm(solve_chol(A1,B1)))
wspp2 = norm(B2 - A2*solve_chol(A2,B2))/(norm(A2)*norm(solve_chol(A2,B2)))
wspp3 = norm(B3 - A3*solve_chol(A3,B3))/(norm(A3)*norm(solve_chol(A3,B3)))
wspp4 = norm(B4 - A4*solve_chol(A4,B4))/(norm(A4)*norm(solve_chol(A4,B4)))
wspp5 = norm(B5 - A5*solve_chol(A5,B5))/(norm(A5)*norm(solve_chol(A5,B5)))
wspp6 = norm(B6 - A6*solve_chol(A6,B6))/(norm(A6)*norm(solve_chol(A6,B6)))
wspp7 = norm(B7 - A7*solve_chol(A7,B7))/(norm(A7)*norm(solve_chol(A7,B7)))
wspp8 = norm(B8 - A8*solve_chol(A8,B8))/(norm(A8)*norm(solve_chol(A8,B8)))

% maksymalny współczynnik wynosi także 10^(-15)



tic; rozwUkladu1(A1, B1); toc
tic; rozwUkladu2(A1,B1); toc
tic; solve_chol(A1,B1); toc
tic; inv(A1)*B1; toc

tic; rozwUkladu1(A2, B2); toc
tic; rozwUkladu2(A2,B2); toc
tic; solve_chol(A2,B2); toc
tic; inv(A2)*B2; toc

tic; rozwUkladu1(A3, B3); toc
tic; rozwUkladu2(A3,B3); toc
tic; solve_chol(A3,B3); toc
tic; inv(A3)*B3; toc

tic; rozwUkladu1(A4, B4); toc
tic; rozwUkladu2(A4,B4); toc
tic; solve_chol(A4,B4); toc
tic; inv(A4)*B4; toc

tic; rozwUkladu1(A5,B5); toc
tic; rozwUkladu2(A5,B5); toc
tic; solve_chol(A5,B5); toc
tic; inv(A5)*B5; toc

tic; rozwUkladu1(A6,B6); toc
tic; rozwUkladu2(A6,B6); toc
tic; solve_chol(A6,B6); toc
tic; inv(A6)*B6; toc

tic; rozwUkladu1(A7,B7); toc
tic; rozwUkladu2(A7,B7); toc
tic; solve_chol(A7,B7); toc
tic; inv(A7)*B7; toc

a=tic; rozwUkladu1(A8,B8); toc
b=tic; rozwUkladu2(A8,B8); toc
c=tic; solve_chol(A8,B8); toc
d=tic; inv(A8)*B8; toc

hold on
plot(1,0.000195,'r:x')
plot(1,0.000979,'g:x')
plot(1,0.000273,'b:x')
plot(1,0.000145,'k:x')

plot(2,0.000167,'r:x')
plot(2,0.000755,'g:x')
plot(2,0.000144,'b:x')
plot(2,0.000086,'k:x')

plot(3,0.000088,'r:x')
plot(3,0.001064,'g:x')
plot(3,0.000257,'b:x')
plot(3,0.000121,'k:x')

plot(4,0.000150,'r:x')
plot(4,0.001196,'g:x')
plot(4,0.000267,'b:x')
plot(4,0.000152,'k:x')

plot(5,0.000195,'r:x')
plot(5,0.004725,'g:x')
plot(5,0.000350,'b:x')
plot(5,0.000136,'k:x')

plot(6,0.000561,'r:x')
plot(6,0.020988,'g:x')
plot(6,0.000855,'b:x')
plot(6,0.000246,'k:x')

plot(7,0.001530,'r:x')
plot(7,0.066086,'g:x')
plot(7,0.002837,'b:x')
plot(7,0.000344,'k:x')

plot(8,0.006781,'r:x')
plot(8,0.313887,'g:x')
plot(8,0.323598,'b:x')
plot(8,0.324577,'k:x')

hold off

ERR1 = abs(solve_chol(A1,B1) - inv(A1)*B1);
imagesc(ERR1);
colorbar;
title('Macierz ERR1');


A =solve_chol(A1,B1) - inv(A1)*B1
B =solve_chol(A2,B2) - inv(A2)*B2
C =solve_chol(A3,B3) - inv(A3)*B3
D =solve_chol(A4,B4) - inv(A4)*B4
E =solve_chol(A5,B5) - inv(A5)*B5
F =solve_chol(A6,B6) - inv(A6)*B6
G =solve_chol(A7,B7) - inv(A7)*B7
H =solve_chol(A8,B8) - inv(A8)*B8

ERR8 = abs(solve_chol(A8,B8) - inv(A8)*B8);
imagesc(ERR8);
colorbar;
title('Macierz ERR8');

