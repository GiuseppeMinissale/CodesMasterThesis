clear

global g b L k c beta1 mu AvHm1 AvHm2 delta beta2
vnorms = [];
vt = [];
energy = []; 
d_energy = []; 

g = 9.81;
L = 1000;
tfin = 200; 
delta = 1;

q = 10;       
N = 2^q;     
dx = L / N;
x = -L / 2 : dx : L / 2 - dx;

b = -3 / 5 + 2 / 5 * sin(2 * pi * x);
H = -b;
AvHm1 = mean(H.^-1); 
AvHm2 = mean(H.^-2);
flucHm1 = (1 ./ H) - AvHm1;
flucflucHm1 = trapz(x, flucHm1); 


c = sqrt(g / AvHm1); 
mu = mean(flucflucHm1.^2)/AvHm1.^2; 
beta2 = delta * AvHm2 / AvHm1;
beta1 = beta2 * c^2;

[eta, q] = initial(x);

Uh = [fft(eta); fft(q)];

CFL = 0.9;

k = -N / 2 : N / 2 - 1;
k = k * 2 * pi / L;
k = fftshift(k);

eta = real(ifft(Uh(1, :)));
q = real(ifft(Uh(2, :)));
h = eta - b;
u = q ./ h;
cmax = max(abs(u) + sqrt(g * h));
dt = CFL * dx / cmax;

t = 0;
n = 0;

while t < tfin
    Uh_new = RK4h(@funchlag, Uh, dt);
    
    [eta, q] = field(Uh_new);
    
    % Calcolo di q_x
    q_x = real(ifft(1i * k .* Uh_new(2, :))); 
    
    E = (c^2 / 2) * eta.^2 + (1 / 2) * q.^2 + (mu / 2) * q_x.^2 + (1 / 6) * beta1 * eta.^3;
    E_total = sum(E) * dx;
    if n > 0
        dE_dt = (E_total - energy(end)) / dt;
        d_energy = [d_energy, dE_dt];
    end
    
    energy = [energy, E_total]; 
    
    norms = cnorms(Uh_new);
    vnorms = [vnorms; norms];
    vt = [vt, t];
    
    if mod(n, 200) == 0
        visualize(x, eta, q, t);
        pause(0.01);
    end

    Uh = Uh_new;
    t = t + dt;
    n = n + 1;
end

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

figure(5)

subplot(3, 1, 1)
plot(vt, energy, 'LineWidth', 1.5)
xlabel('Time')
ylabel('E(t)')
title('E(t) (CFL=0.9)')
grid on

subplot(3, 1, 2)
plot(vt, energy, 'LineWidth', 1.5)
xlabel('Time')
ylabel('E(t)')
title('E(t) on a bigger scale (CFL=0.9)')
grid on

ylim([min(energy) - 10^-2, max(energy) + 10^-2])  % Moltiplica i limiti dell'asse y per 10^2

subplot(3, 1, 3)
plot(vt(2:end), d_energy, 'LineWidth', 1.5)
xlabel('Tempo')
ylabel('dE/dt')
title('Derivative of E (CFL=0.9)')
grid on

% --- FUNZIONI AUSILIARIE ---
function u = RK4h(funch, u, dt)
    K1 = dt * funch(u);
    K2 = dt * funch(u + K1 / 2);
    K3 = dt * funch(u + K2 / 2);
    K4 = dt * funch(u + K3);
    u = u + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
end

function fh = funchlag(Uh)
    global g b k c beta1 beta2 mu
    etah = Uh(1, :);
    qh = Uh(2, :);
    eta = ifft(etah);
    q = ifft(qh);
    
    eta_x = ifft(1i * k .* etah);
    q_x = ifft(1i * k .* qh);
    
    fh = zeros(size(Uh));
    fh(1, :) = -1i * k .* qh; % Derivata spaziale di q
    fh(2, :) = (-c^2 * 1i *k.*etah - 2 * beta2 * fft(q .* q_x) - beta1 * fft(eta .* eta_x)) ./ (1 + mu * k.^2);
end

function [eta, q] = initial(x)
    global L
    eta = 0.025 * exp(-x.^2 / 9);
    q = zeros(size(x));
end

function visualize(x, eta, q, t)
    global b
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