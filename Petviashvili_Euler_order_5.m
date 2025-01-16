% Parametri globali
global g ps vs m mu Km1 V cs2 delta

% Costanti di base
g = 1.4;
ps = 1;
vs = 1;
m = 1 + 1/g;
delta = 1;
mu = 0.01171875;
Km1 = 1;   % <K^{-1}>
V = 1.2;
cs2 = g * ps / vs;
zeta = 0.0004577636718750001;
kappa = 1.4;
nu = delta^4 * (zeta / Km1^3 - mu^2);

% Definire intervallo spaziale per calcolo
x = linspace(-0.2, 0.2, 80);
y = linspace(-0.1, 0.1, 80);

% Funzioni per le quantità calcolate
G = @(p) cs2 * (p / ps).^m;
G1 = @(p) m ./ p .* G(p);
GG = @(p) p / (m + 1) .* G(p);

Ft = @(x, y) m / ps * G(ps) / (2 * V * Km1) * y.^2 ...
    - 1 / (mu * delta^2 * V^2) * (((ps + V * x) .* G(ps + V * x) - ps * G(ps)) / (Km1 * (m + 1) * V) - V^2 * x);

F = @(x, y) m * cs2 / (2 * ps * V) * y.^2 ...
    - 1 / (mu * delta * V^2) * (cs2 * ((1 + V * x / ps).^m .* (ps + V * x) - ps) / ((m + 1) * Km1 * V) - V^2 * x);

df = F(x, y) - Ft(x, y);

Ftt = @(x, y) G1(ps) / (2 * Km1 * V) * y.^2 ...
    - G1(ps) / (2 * delta^2 * mu * V * Km1 * m) * x.^2;

[X, Y] = meshgrid(x, y);
Fx = Y;
Fy = Ftt(X, Y);
Fx = Fx ./ sqrt(Fx.^2 + Fy.^2);
Fy = Fy ./ sqrt(Fx.^2 + Fy.^2);

quiver(X, Y, Fx, Fy)

% Calcolo della matrice per il separatore di fase
beta = (1 - G(ps) / V^2 / Km1) / (delta^2 * mu);
A = [0 1; beta 0];
[vec, lam] = eig(A);

% Valori iniziali
epsilon = 0.06e-6;
u0 = [epsilon * vec(1, 1); epsilon * vec(2, 1)];
tmax = 22;

options = odeset('reltol', 1e-16, 'abstol', 1e-14, 'InitialStep', 1e-14, 'MaxStep', 1e-3);
T_nodes = linspace(0, tmax ,10000);
% Risolvi il sistema ODE rispetto a xi
[T, U] = ode78('udot3', T_nodes, u0, options);

% Plot della soluzione nella fase
figure(5)
subplot(1, 2, 1)
plot(U(:, 1), U(:, 2), u0(1), u0(2), 'o', 'linewidth', 3)
xlabel('$u$', 'fontsize', 32, 'interpreter', 'latex')
ylabel('$u''$', 'fontsize', 32, 'interpreter', 'latex')
fontsize(24, "points")

subplot(1, 2, 2)
plot(T, U(:, 1), 'linewidth', 3)
axis([0 20 0 0.1])
xlabel('$\xi$', 'fontsize', 32, 'interpreter', 'latex')
ylabel('$u(\xi)$', 'fontsize', 32, 'interpreter', 'latex')
fontsize(24, "points")

% Parametri per la soluzione Fourier
k1 = G1(ps)/(Km1*m*V*nu*2);
c1 = (1/nu-G(ps)/(Km1*V^2*nu));
k2 = -delta^2*mu*G1(ps)/(2*Km1*V*nu);
c2 = -(delta^2*mu/nu); 
max_iter = 1000;  % Numero massimo di iterazioni
  % Numero di punti nella soluzione iniziale U

% Funzioni per la costruzione della matrice e calcolo di T
function L = costruisci_matrice_L(N, T, c1, c2)
    I = speye(N);
    tmin=min(T);
    tmax=max(T);
    l_T=abs(tmax-tmin);
    factor2 = (N-1)^2/l_T^2;
    factor4= (N-1)^4/l_T^4;
    
    
% Creazione della matrice a bande usando spdiags
    D2 = spdiags([1,-2,1], [-1, 0, 1], N, N);


    D2(N,1) = 1;
    D2(1,N) = 1;
    D4=matrice_derivata_quarta(N);
    
   

   
    L = c1 * I + c2 * factor2 * D2+ factor4* D4;
   
end


function D4 = matrice_derivata_quarta(N)
   
    
    % Numero di punti
    
    
    % Costruzione della matrice di differenze finite centrali per la derivata quarta
    D4 = spdiags([1,-4,6,-4,1], [-2,-1, 0, 1,2], N, N);
    
    D4(2,1)=-4;
    D4(2,3)=-4;
    D4(2,4)=1;
    D4(N-1,N)=-4;
    D4(N-1,N-2)=-4;
    D4(N-1,N-3)=1;
    D4(1,2)=-4;
    D4(1,3)=1;
    D4(1, N-1) = 1;    % Collegamento per x_0-2 (periodicità)
    D4(1, N) = -4;     % Collegamento per x_0-1 (periodicità)
    D4(2, N) = 1;      % Collegamento per x_1-2 (periodicità)
    D4(N,N-2)=1;
    D4(N,N-1)=-4;
    D4(N, 1) = -4;     % Collegamento per x_N+1 (periodicità)
    D4(N, 2) = 1;      % Collegamento per x_N+2 (periodicità)
    D4(N-1, 1) = 1;    % Collegamento per x_N-1+2 (periodicità) 
    
 
end


function O_u = calcola_O(u, T, k1, k2)
    N = length(u);
    O_u = zeros(N, 1);
    tmin=min(T);
    tmax=max(T);
    l_T=abs(tmax-tmin);
    factor = N^2/(4*l_T^2);

    for i = 1:N
        if i == 1
            O_u(i) = k1 * u(i)^2 + k2 * factor * (u(i + 1) - u(N))^2;
        elseif i == N
            O_u(i) = k1 * u(i)^2 + k2 * factor * (u(1)- u(i - 1))^2;
        else
            O_u(i) = k1 * u(i)^2 + k2 * factor * (u(i + 1) - u(i - 1))^2;
        end
    end
end




function u_final = iterative_solution_Petviashvili(U, T, c1, c2, k1, k2, max_iter, tol)
    N=length(U);
    L = costruisci_matrice_L(N, T, c1, c2);
    L_inv = inv(L);
    u_n = U(:, 1);  % Condizione iniziale
    figure;
    hold on;
    xlabel('\xi');
    ylabel('u(\xi)');
    title('u determinata tramite il metodo di Petviashvili');
    plot_handle = plot(T, u_n);
    pause(1);
    for n = 1:max_iter
        O_un = calcola_O(u_n, T,  k1, k2);
        L_un = L * u_n;
        s_n = (norm(L_un, 2) / norm(O_un, 2))^2;
        u_nplus1 = L\(s_n * O_un);
        diff_norm = norm(u_nplus1 - u_n, inf);
        fprintf('Iterazione %d: norma della differenza = %.5f\n', n, diff_norm);
        residual_norm = norm(L_un - O_un, inf);
        fprintf('Iterazione %d: norma di L*u_n - T(u_n) = %.5f\n', n, residual_norm);
        set(plot_handle, 'YData', u_nplus1);
        drawnow;
        pause(0.1);
        u_n = u_nplus1;
        % Criterio di arresto basato sulla differenza tra u_n+1 e u_n
        if diff_norm < tol
            fprintf('Convergenza raggiunta in %d iterazioni\n', n);
            u_final = u_nplus1;
            return;
        end
        
        % Aggiorna u_n per la prossima iterazione
        u_n = u_nplus1;
    end
    u_final = u_n;
end
tol=0.000001;
u_final = iterative_solution_Petviashvili(U.^2, T, c1, c2, k1, k2, max_iter, tol);

figure
hold on
plot(T, u_final(:,1), 'b', 'LineWidth', 1)
title('Travelling wave obtained with Petviashvili method')
hold off



figure(2)
hold on
plot(T, U(:,1), 'r', 'LineWidth', 0.1) % Onda viaggiante Ftt1 (rosso)
plot(T, u_final(:,1), 'b', 'LineWidth', 0.1) % Onda viaggiante Ftt2 (blu)
xlabel('$\xi$', 'fontsize', 32, 'interpreter', 'latex')
ylabel('$u(\xi)$', 'fontsize', 32, 'interpreter', 'latex')
legend('U','u')
title('Comparison','U obtained with classic approach vs u obtained with Petviashvili method','FontWeight','bold')
hold off


figure
hold on
semilogy(T,abs(u_final(:,1)-U(:,1)), 'r', 'LineWidth', 1) 
xlabel('$\xi$', 'fontsize', 32, 'interpreter', 'latex')
ylabel('$|U-u|$', 'fontsize', 32, 'interpreter', 'latex')
legend('|U-u|')
title('|U-u|','Difference between the two travelling waves', 'FontWeight','bold')
hold off
