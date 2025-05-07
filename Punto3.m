clc
clear all
close all

% Intervallo di frequenze (omega) desiderato, con passo 0.01
omega = 0:0.001:100;

% Inizializzare i vettori per il modulo e la fase per entrambe le funzioni
magnitude_real = zeros(size(omega));
phase_real = zeros(size(omega));
magnitude_approx = zeros(size(omega));
phase_approx = zeros(size(omega));
magnitude_Lreal = zeros(size(omega));
phase_Lreal = zeros(size(omega));
magnitude_Lapprox = zeros(size(omega));
phase_Lapprox = zeros(size(omega));

% Calcolare il modulo e la fase per ciascuna frequenza per entrambe le funzioni
for k = 1:length(omega)
    s = 1i * omega(k);
    kp = 0.006646622313;
    omen = 2*pi*50;
    An = 0.29337;
    kstatic = omen/An;
    T = 0.3507;
    TA = 12;
    Beta = 0.3811;
    mu=omen/(340e6);

    % Calcolo per Greale
    H_real = (1 - 2* Beta * tanh(0.91 * s)) / (1 + Beta * tanh(0.91 * s));
    magnitude_Hreal(k) = abs(H_real);
    phase_Hreal(k) = angle(H_real) * (180 / pi);
    
    % Calcolo per Gapprox
    H_approx = (1 - 2*T * s) / (1 + T * s);
    magnitude_Happrox(k) = abs(H_approx);
    phase_Happrox(k) = angle(H_approx) * (180 / pi);
    
    % Calcolo per Lreale
    factor = kp * kstatic * ((20*s+1)/(20*s)) * (1/(TA*s+1)) * (1/(0.2*s+1));
    L_real = H_real * factor;
    magnitude_Lreal(k) = abs(L_real);
    phase_Lreal(k) = angle(L_real) * (180 / pi); % Convertire la fase in gradi
    
    % Calcolo per Lapprox
    L_approx = H_approx * factor;
    magnitude_Lapprox(k) = abs(L_approx);
    phase_Lapprox(k) = angle(L_approx) * (180 / pi); % Convertire la fase in gradi

     % Calcolo Ssapprox
    Ss_approx = mu*(1/(TA*s+1))*(1/(1+L_approx));
    magnitude_Ssapprox(k) = abs(Ss_approx);
    phase_Ssapprox(k) = angle(Ss_approx) * (180 / pi);

    % Calcolo Ssreale
    Ss_reale = mu*(1/(TA*s+1))*(1/(1+L_real));
    magnitude_Ssreale(k) = abs(Ss_reale);
    phase_Ssreale(k) = angle(Ss_reale) * (180 / pi);
end


% Plottare il diagramma di modulo (in dB) per Lreale e Lapprox
figure (1);
subplot(2, 1, 1);
semilogx(omega, 20 * log10(magnitude_Lreal), 'b', omega, 20 * log10(magnitude_Lapprox), 'r--');
grid on;
title('Diagramma di Bode - Modulo di Lreale e Lapprox');
xlabel('Frequenza (rad/s)');
ylabel('Ampiezza (dB)');
legend('Lreale', 'Lapprox');

hold
% Plottare il diagramma di fase (in gradi) per Lreale e Lapprox
subplot(2, 1, 2);
semilogx(omega, phase_Lreal, 'b', omega, phase_Lapprox, 'r--');
grid on;
title('Diagramma di Bode - Fase di Lreale e Lapprox');
xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
legend('Lreale', 'Lapprox');

hold
% Plottare il diagramma di modulo (in dB) per Ssreale e Ssapprox
figure (2);
subplot(2, 1, 1);
semilogx(omega, 20 * log10(magnitude_Ssreale), 'b', omega, 20 * log10(magnitude_Ssapprox), 'r--');
grid on;
title('Diagramma di Bode - Modulo di Ssreale e Ssapprox');
xlabel('Frequenza (rad/s)');
ylabel('Ampiezza (dB)');
legend('Ssreale', 'Ssapprox');

hold
% Plottare il diagramma di fase (in gradi) per Lreale e Lapprox
subplot(2, 1, 2);
semilogx(omega, phase_Ssreale, 'b', omega, phase_Ssapprox, 'r--');
grid on;
title('Diagramma di Bode - Fase di Ssreale e Ssapprox');
xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
legend('Ssreale', 'Ssapprox');

hold

% Plottare il diagramma di modulo (in dB) per Hreale e Happrox
figure (3);
subplot(2, 1, 1);
semilogx(omega, 20 * log10(magnitude_Hreal), 'b', omega, 20 * log10(magnitude_Happrox), 'r--');
grid on;
title('Diagramma di Bode - Modulo di Hreale e Happrox');
xlabel('Frequenza (rad/s)');
ylabel('Ampiezza (dB)');
legend('Hreale', 'Happrox');

hold
% Plottare il diagramma di fase (in gradi) per Hreale e Happrox
subplot(2, 1, 2);
semilogx(omega, phase_Hreal, 'b', omega, phase_Happrox, 'r--');
grid on;
title('Diagramma di Bode - Fase di Hreale e Happrox');
xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
legend('Hreale', 'Happrox');

hold
%% Point 7

% Intervallo di frequenze (omega) desiderato, con passo 0.01
omega = 0:0.001:100;

% Inizializzare i vettori per il modulo e la fase per entrambe le funzioni
magnitude_real = zeros(size(omega));
phase_real = zeros(size(omega));
magnitude_approx = zeros(size(omega));
phase_approx = zeros(size(omega));
magnitude_Lreal = zeros(size(omega));
phase_Lreal = zeros(size(omega));
magnitude_Lapprox = zeros(size(omega));
phase_Lapprox = zeros(size(omega));

% Calcolare il modulo e la fase per ciascuna frequenza per entrambe le funzioni
for k = 1:length(omega)
    s = 1i * omega(k);
    kp = 0.006646622313;
    omen = 2*pi*50;
    An = 0.29337/5;
    kstatic = omen/An;
    T = 0.3507/5;
    TA = 12*5;
    Beta = 0.3861/5;
    alpha = 1*5;
    mu=omen/(340e6)*5;

    % Calcolo per Greale
    H_real = (1 - 0.2*Beta * tanh(0.91 * s)) / (1 + Beta* tanh(0.91 * s));
    
    % Calcolo per Gapprox
    H_approx = (1 - 2*T * s) / (1 + T * s);
    
    % Calcolo per Lreale
    factor = kp * kstatic * ((20*s+1)/(20*s)) * (1/(TA*s+alpha)) * (1/(0.2*s+1));
    L_real = H_real * factor;
    magnitude_Lreal(k) = abs(L_real);
    phase_Lreal(k) = angle(L_real) * (180 / pi); % Convertire la fase in gradi
    
    % Calcolo per Lapprox
    L_approx = H_approx * factor;
    magnitude_Lapprox(k) = abs(L_approx);
    phase_Lapprox(k) = angle(L_approx) * (180 / pi); % Convertire la fase in gradi

    % Calcolo Ssapprox
    Ss_approx = mu*(1/(TA*s+alpha))*(1/(1+L_approx));
    magnitude_Ssapprox(k) = abs(Ss_approx);
    phase_Ssapprox(k) = angle(Ss_approx) * (180 / pi);

    % Calcolo Ssreale
    Ss_reale = mu*(1/(TA*s+alpha))*(1/(1+L_real));
    magnitude_Ssreale(k) = abs(Ss_reale);
    phase_Ssreale(k) = angle(Ss_reale) * (180 / pi);
end

% Plottare il diagramma di modulo (in dB) per Lreale e Lapprox
figure (1);
subplot(2, 1, 1);
semilogx(omega, 20 * log10(magnitude_Lreal), omega, 20 * log10(magnitude_Lapprox), '--');
grid on;
title('Diagramma di Bode - Modulo di Lreale e Lapprox');
xlabel('Frequenza (rad/s)');
ylabel('Ampiezza (dB)');
legend('Lreale*', 'Lapprox*');

% Plottare il diagramma di fase (in gradi) per Lreale e Lapprox
subplot(2, 1, 2);
semilogx(omega, phase_Lreal, omega, phase_Lapprox, '--');
grid on;
title('Diagramma di Bode - Fase di Lreale e Lapprox');
xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
legend('Lreale*', 'Lapprox*');


% Plottare il diagramma di modulo (in dB) per Ssreale e Ssapprox
figure (2);
subplot(2, 1, 1);
semilogx(omega, 20 * log10(magnitude_Ssreale), omega, 20 * log10(magnitude_Ssapprox), '--');
grid on;
title('Diagramma di Bode - Modulo di Lreale e Lapprox');
xlabel('Frequenza (rad/s)');
ylabel('Ampiezza (dB)');
legend('Ssreale*', 'Ssapprox*');

% Plottare il diagramma di fase (in gradi) per Lreale e Lapprox
subplot(2, 1, 2);
semilogx(omega, phase_Ssreale, omega, phase_Ssapprox, '--');
grid on;
title('Diagramma di Bode - Fase di Lreale e Lapprox');
xlabel('Frequenza (rad/s)');
ylabel('Fase (gradi)');
legend('Ssreale*', 'Ssapprox*');