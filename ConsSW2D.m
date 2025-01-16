clear

global g b L k mu a1 a2 delta AvH tmu
vnorms = [];
vt = [];
energy = []; % Per memorizzare i valori di E(t)
energy_E0=[];
energy_E1=[];
energy_E3=[];
energy_E4=[];
d_energy = []; % Per memorizzare la derivata di E(t)

% Parametri globali e fisici
g = 9.81; % Accelerazione dovuta alla gravità
L = 1000; % Lunghezza del dominio spaziale
tfin = 200; % Tempo finale
delta = 1;

% Parametri spaziali e temporali
q = 10;       % Potenza di due (dimensioni gestibili)
N = 2^q;      % Numero di intervalli
dx = L / N;
x = -L / 2 : dx : L / 2 - dx;

% Profondità del fondale periodico
b = -1+ 3 / 10 * sin(2 * pi * x);
H = -b;
AvH = mean(H); % Approssimazione numerica dell'integral
a1=g*AvH;
a2=delta/AvH;
mu=0.0483; %calcoli fatti con python a partire dal codice del prof. Ketcheson trovato online
tmu=delta^2*mu./AvH;
disp(tmu);




[eta, q] = initial(x);

% FFT della condizione iniziale
Uh = [fft(eta); fft(q)];

% Parametri numerici
CFL = 0.9; % Numero di Courant-Friedrichs-Lewy

% Spazio dei numeri d'onda
k = -N / 2 : N / 2 - 1;
k = k * 2 * pi / L;
k = fftshift(k);

% Calcolo passo temporale iniziale
eta = real(ifft(Uh(1, :)));
q = real(ifft(Uh(2, :)));
h = eta - b;
u = q ./ h;
cmax = max(abs(u) + sqrt(g * h));
disp(cmax);
dt = CFL * dx / cmax;

% Ciclo temporale
t = 0;
n = 0;

while t < tfin
    % Runge-Kutta di quarto ordine
    Uh_new = RK4h(@funchlag, Uh, dt);
    
    % Recupero eta e q
    [eta, q] = field(Uh_new);
    
    % Calcolo di q_x
    q_x = real(ifft(1i * k .* Uh_new(2, :))); % Derivata spaziale di q

    % Calcolo di E(t)
    E0 =1/2*q.^2+a1/2*eta.^2+tmu/2*q_x.^2; %E fino a delta^0
    E1 =1/2*q.^2+a1/2*eta.^2+tmu/2*q_x.^2-a1*a2/6*eta.^3; %E fino a delta
    E2 =1/2*q.^2+a1/2*eta.^2+tmu/2*q_x.^2-a1*a2/6*eta.^3+(a1*a2^2)/12*eta.^4; %E fino a delta^2
    E3 =1/2*q.^2+a1/2*eta.^2+tmu/2*q_x.^2-a1*a2/6*eta.^3+(a1*a2^2)/12*eta.^4-(a1*a2^3)/20*eta.^5; %E fino a delta^3
    E4 =1/2*q.^2+a1/2*eta.^2+tmu/2*q_x.^2-a1*a2/6*eta.^3+(a1*a2^2)/12*eta.^4-(a1*a2^3)/20*eta.^5+(a1*a2^4)/30*eta.^6; %E fino a delta^4
    E_total= sum(E2)*dx;
    E_total0= sum(E0)*dx;
    E_total1= sum(E1)*dx;
    E_total3= sum(E3)*dx;
    E_total4 = sum(E4) * dx; 
    
    
    % Calcolo della derivata temporale di E(t)
    if n > 0
        dE_dt = (E_total4 - energy_E4(end)) / dt;
        d_energy = [d_energy, dE_dt];
    end
    
    energy = [energy, E_total]; 
    energy_E4 = [energy_E4, E_total4];
    energy_E3 = [energy_E3, E_total3];
    energy_E1 = [energy_E1, E_total1];
    energy_E0 = [energy_E0, E_total0];


    % Calcolo delle norme
    norms = cnorms(Uh_new);
    vnorms = [vnorms; norms];
    vt = [vt, t];
    
    % Visualizzazione
    if mod(n, 200) == 0
        visualize(x, eta, q, t);
        pause(0.01);
    end
    
    % Aggiornamento
    Uh = Uh_new;
    t = t + dt;
    n = n + 1;

end

% Plot delle norme
figure(4)
subplot(2, 1, 1)
plot(vt, vnorms(:, 1), vt, vnorms(:, 2))
xlabel('Tempo')
ylabel('Derivate di \eta')
title('Derivate di \eta')

subplot(2, 1, 2)
plot(vt, vnorms(:, 3), vt, vnorms(:, 4))
xlabel('Tempo')
ylabel('Derivate di q')
title('Derivate di q')

% Plot dell'energia e della sua derivata
figure(5)

% Energia totale su scala naturale
subplot(3, 1, 1)
plot(vt, energy_E4, 'LineWidth', 1.5)
xlabel('Time')
ylabel('E4(t)')
title('E4(t) (CFL=0.45)')
grid on

% Energia totale con cambiamento della scala (cambia solo l'asse y)
subplot(3, 1, 2)
plot(vt, energy_E4, 'LineWidth', 1.5)
xlabel('Time')
ylabel('E4(t)')
title('E4(t) on a bigger scale (CFL=0.45)')
grid on

% Modifica l'asse Y per rappresentarlo su un ordine di grandezza maggiore
ylim([min(energy) - 2*10^-2, max(energy) + 2*10^-2])  % Moltiplica i limiti dell'asse y per 10^2

% Derivata dell'energia totale
subplot(3, 1, 3)
plot(vt(2:end), d_energy, 'LineWidth', 1.5)
xlabel('Time')
ylabel('dE4/dt')
title('Derivative of E4 (CFL=0.45)')
grid on

% Calcolo delle differenze relative finali
diff_E0 = (energy_E0(1) - energy_E0(end)) / energy_E0(end);
diff_E1 = (energy_E1(1) - energy_E1(end)) / energy_E1(end);
diff_E2 = (energy(1) - energy(end)) / energy(end);
diff_E3 = (energy_E3(1) - energy_E3(end)) / energy_E3(end);
diff_E4 = (energy_E4(1) - energy_E4(end)) / energy_E4(end);

% Stampa delle differenze relative
fprintf('Relative difference for E_0: %.6f\n', diff_E0);
fprintf('Relative difference for E_1: %.6f\n', diff_E1);
fprintf('Relative difference for E_2: %.6f\n', diff_E2);
fprintf('Relative difference for E_3: %.6f\n', diff_E3);
fprintf('Relative difference for E_4: %.6f\n', diff_E4);



% --- FUNZIONI AUSILIARIE ---
function u = RK4h(funch, u, dt)
    K1 = dt * funch(u);
    K2 = dt * funch(u + K1 / 2);
    K3 = dt * funch(u + K2 / 2);
    K4 = dt * funch(u + K3);
    u = u + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
end

function fh = funchlag(Uh)
    global k a1 a2 tmu
    etah = Uh(1, :);
    qh = Uh(2, :);
    eta = ifft(etah);
    q = ifft(qh);
    
    % Derivate spaziali
    eta_x = ifft(1i * k .* etah);
    q_x = ifft(1i * k .* qh);
    
    % Sistema delle equazioni
    fh = zeros(size(Uh));
    fh(1, :) = -1i*k.*qh-a2*fft(eta_x.*q)-a2*fft(eta.*q_x); 
    fh(2, :) = (-a1*1i*k.*etah-a2*fft(q.*q_x)) ./ (1 + tmu*k.^2);
end

function [eta, q] = initial(x)
    eta = 0.05 * exp(-x.^2 / 25);
    q = zeros(size(x));
end

function visualize(x, eta, q, t)
    subplot(2, 1, 1)
    plot(x, eta, 'b', 'LineWidth', 1.5)
    xlabel('x')
    ylabel('\eta')
    title(['Elevazione della superficie dell''acqua a t = ', num2str(t)])
    
    subplot(2, 1, 2)
    plot(x, q, 'r', 'LineWidth', 1.5)
    xlabel('x')
    ylabel('q')
    title(['Flusso di quantità di moto a t = ', num2str(t)])
end

function [eta, q] = field(Uh)
    eta = real(ifft(Uh(1, :)));
    q = real(ifft(Uh(2, :)));
end

function norms = cnorms(Uh)
    global k
    qh = Uh(1, :);
    q1 = max(real(ifft(1i * k .* qh)));
    q2 = max(real(ifft((1i * k).^2 .* qh)));
    etah = Uh(2, :);
    eta1 = max(real(ifft(1i * k .* etah)));
    eta2 = max(real(ifft((1i * k).^2 .* etah)));
    norms = [q1, q2, eta1, eta2];
end
