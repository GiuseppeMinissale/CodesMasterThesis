% Parametri
c = 1;         % Velocità dell'onda
N = 200;       % Numero di punti
Lx = 20;       % Lunghezza del dominio
x = linspace(-Lx/2, Lx/2, N)';
dx = x(2) - x(1);

% Operatore differenziale (Matrice a bande)
Dxx = spdiags([ones(N, 1), -2 * ones(N, 1), ones(N, 1)], -1:1, N, N) / dx^2;
Dxx(1, N) = 1 / dx^2;   % Condizioni periodiche
Dxx(N, 1) = 1 / dx^2;

L = -c * speye(N) + Dxx; % Operatore lineare
L_inv = inv(L);          % Inversione del lineare

% Condizione iniziale (ad esempio un impulso gaussiano)
u_n = exp(-x.^2);

% Parametri del metodo di Petviashvili
k1 = -3;  % Coefficiente del termine non lineare
max_iter = 1000;
tol = 1e-6;

% Iterazione di Petviashvili
for n = 1:max_iter
    % Calcola il termine non lineare
    T_un = k1 * u_n.^2;
    
    % Calcola s_n
    L_un = L * u_n;
    s_n = (norm(L_un, 2) / norm(T_un, 2))^2;
    
    % Aggiorna u_n+1
    u_nplus1 = L_inv * (s_n * T_un);
    
    % Criterio di arresto
    diff_norm = norm(u_nplus1 - u_n, inf);
    fprintf('Iterazione %d: norma della differenza = %.5e\n', n, diff_norm);
    if diff_norm < tol
        fprintf('Convergenza raggiunta in %d iterazioni\n', n);
        break;
    end
    
    u_n = u_nplus1; % Aggiorna
    
end

% Risultati
figure;
plot(x, u_n, 'LineWidth', 2);
title('Travelling wave computed using Petviashvili method');
xlabel('\xi');
ylabel('u(\xi)');
grid on;
