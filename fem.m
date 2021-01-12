% liczba elementów
N = 20;

[M,func] = fill_matrix(N); 
u = get_u(M,func,N); 

% disp(M);

% wykres funkcji u(x)
fplot(u, [0,2]);
ylim([0,30]);
title('Przybliżenie u(x)');
figure

% wykres funkcji bazowych (kształtu)
fplot(func{1} ,[0 2])
title('Funkcje bazowe');
hold on
for k = 2:N+1
    fplot(func{k}, [0 2])
end
hold off


% f. zwracająca przybliżenie funkcji u, będącej kombinacją liniową f. bazowych
function u = get_u(M,func,N)
    B = zeros(N,1);
    B(1) = 30;
    X = linsolve(M,B);
    
    u = 0;
    for k = 1:N
        u = u + func{k}.*X(k);
    end  
end
        

% f. zwracająca odpowiednią macierz oraz funkcje bazowe
function [M,func] = fill_matrix(N)
    syms E(x)
    E(x) = piecewise(and(0 <= x, x <= 1), 3, and(1 < x, x <= 2), 5, 0);
    
    seg = 2/N;  % długość jednego przedziału
    
    M = zeros(N,N);
    for i = 1:N
        for j = i:N
            if i == 1
                u(x) = piecewise(and(0 <= x, x <= seg), (seg - x)/seg, 0);

                % M(1,1)
                if j == 1
                    func{j} = u(x);      % cell array
                    M(i,j) = E(0)*u(0)*u(0) - gaussian_quad(E*diff(u)*diff(u),0,seg);
                    
                % 1 wiersz
                else
                    v(x) = piecewise(and((j-2)*seg <= x, x <= (j-1)*seg), (x - (j-2)*seg)/seg, and((j-1)*seg < x, x <= j*seg), (j*seg - x)/seg, 0);
                    M(i,j) = E(0)*u(0)*v(0) - gaussian_quad(E*diff(u)*diff(v),0,seg);
                    M(j,i) = M(i,j);
                    
                    func{j} = v(x);    
                end
           
            else 
                u(x) = func{j};
                v(x) = func{i};

                if i == j
                    M(i,j) = E(0)*u(0)*v(0) - gaussian_quad(E*diff(u)*diff(v),(i-2)*seg,i*seg);
                elseif j - i == 1
                    M(i,j) = E(0)*u(0)*v(0) - gaussian_quad(E*diff(u)*diff(v),(i-1)*seg,i*seg);
                else
                   M(i,j) = 0;
                end
                
                M(j,i) = M(i,j);
            end     
        end
    end
    
    % ostatnia funkcja
    v(x) = piecewise(and((2 - seg) <= x, x <= 2), (x - (2 - seg))/seg, 0);
    func{N+1} = v(x);
end


% int(f, a, b)
% https://en.wikipedia.org/wiki/Gaussian_quadrature
function q = gaussian_quad(f,a,b)

%     points2 = [-1/sqrt(3) 1/sqrt(3)]; % dwa punkty kwadratury
%     weights2 = [1 1];  % wagi

    points4 = [sqrt(3/7-2/7*sqrt(6/5))  -sqrt(3/7-2/7*sqrt(6/5))  sqrt(3/7+2/7*sqrt(6/5))  -sqrt(3/7+2/7*sqrt(6/5))];
    weights4 = [(18+sqrt(30))/36  (18+sqrt(30))/36  (18-sqrt(30))/36  (18-sqrt(30))/36];
    
    % iloczyn skalarny wektorów da szukaną sumę weight*point
    % przesunięcie punktów kwadratury wynika z konieczności
    % przyjęcia granic całkowania jako [-1 1]
    q = ((b-a)/2).*dot(f(points4.*(b-a)/2 + (a+b)/2),weights4);
end
