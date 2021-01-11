% coś tu kurw mocno nie gra : \


% liczba elementów
N = 3;

% rozwiązanie równania liniowego (macierz kwadratowa)
A = [-3/2  9/2   0;
      9/2 -21/2  6;
       0    6   -27/2];
B = [30; 0; 0];
X = linsolve(A,B);


syms E(x)
E(x) = piecewise(and(0 <= x, x <= 1), 3, and(1 < x, x <= 2), 5, 0);

syms e0(x)
e0(x) = piecewise(and(x >= 0, x <= 2/3), -3/2.*x + 1, 0);
e0_diff = diff(e0);

syms e1(x)
e1(x) = piecewise(and(x >= 0, x <= 2/3), 3/2*x, and(x > 2/3, x <= 4/3), -3/2*x+2, 0);
e1_diff = diff(e1);

syms e2(x)
e2(x) = piecewise(and(x >= 0, x <= 2/3), 0, and(x > 2/3, x <= 4/3), 3/2*x-1, -3/2*x + 3);
e2_diff = diff(e2);

syms e3(x)
e3(x) = piecewise(and(x >= 0, x <= 4/3), 0, 3/2*x - 2); 
e3_diff = diff(e3);

% TODO
M = zeros(3,3);
for i = 1:N
    for j = 1:N
    M(i,j) = i;
    end
end
disp(M)


A = [-3/2  9/2   0;
      9/2 -21/2  6;
       0    6   -27/2];
 
B = [30; 0; 0];
X = linsolve(A,B);



fplot(@(x) e0(x)*X(1) + e1(x)*X(2) + e2(x)*X(3) ,[0 2])
ylim([0,30])

figure

% wykres funkcji bazowych (kształtu)
fplot(@(x) e0(x) ,[0 2])
hold on
fplot(@(x) e1(x) ,[0 2])
fplot(@(x) e2(x) ,[0 2])
fplot(@(x) e3(x) ,[0 2])
hold off

disp(double(gaussian_quad(E.*e0_diff.*e0_diff,0,2/3)))




% int(f, a, b)
% https://en.wikipedia.org/wiki/Gaussian_quadrature
function q = gaussian_quad(f,a,b)
       
    points2 = [-1/sqrt(3) 1/sqrt(3)]; % dwa punkty kwadratury
    weights2 = [1 1];  % wagi
    % iloczyn skalarny wektorów da szukaną sumę weight*point
    % przesunięcie punktów kwadratury wynika z konieczności
    % przyjęcia granic całkowania jako [-1 1]
    q = ((b-a)/2).*dot(f(points2.*(b-a)/2 + (a+b)/2),weights2);
end
