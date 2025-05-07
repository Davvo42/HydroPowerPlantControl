% PUNTO 1

% Calcola antitrasformata
c = 1321;
syms s t 
K_s = -4002.5/s*(1/(0.614*exp(-0.908*s) + 1.386*exp(0.908*s)));  % Definisci la funzione K(s)
K_t = ilaplace(K_s, s, t);  % Calcola la trasformata inversa
disp(K_t);  % Mostra il risultato


% f_t_numeric = matlabFunction(dw);
% t_values = 0:0.01:15;
% f_values = f_t_numeric(t_values);
% % Plot 
% plot(t_values, f_values);
% xlabel('Tempo t');
% ylabel('f(t)');
% title('Grafico della funzione f(t)');
% grid on;

% PUNTO 2: IMPLEMENTAZIONE PI

% %parameters
% beta = 0.2922;
% alfa = 1;
% Ti = 20;
% wc = 0.595;
% kp = 6.82;
% 
% % 0D model
% s = tf('s');
% Ga  = tf(1, [0.2 1]);
% Ghyd = tf([-0.7 1],[0.35 1]);
% G3 = tf(1, [12 1]);
% G = Ga*Ghyd*G3;
% R = tf([kp*Ti kp], [Ti 0]);
% L = G*R;
% [Gain_margin,Phase_margin,Wcg, W_cross] = margin(L);
% F = L/(1+L);
% 
% % 1D model
% omega = 0:0.01:100;
% 
% %G1
% g1 = 1./(0.2*1i*omega +1);
% 
% 
% 
% %G2
% b = 0.3861;
% a = b*2;
% c = 0.91;
% y = (1 - a * tanh(c *1i*omega)) ./ (1 + b * tanh(c *1i*omega));
% % Calcola il modulo della funzione di trasferimento
% mod_y = abs(y);
% ph_y = angle(y);
% 
% % Plot del modulo in funzione di omega
% figure;
% plot(20*log10(omega), 20*log10(mod_y)); % Conversione in decibel
% xlabel('Frequency (rad/s)');
% ylabel('Magnitude (dB)');
% title('Bode Plot of the Magnitude G2');
% grid on;
% 
% % figure;
% % plot(20*log10(omega), 180/pi*ph_y); % Conversione in decibel
% % xlabel('Frequency (rad/s)');
% % ylabel('Phase');
% % title('Bode Plot of the Phase G2');
% % grid on;
% % 
% % %comparison
% % g2 = (-0.7*1i*omega+1)./(0.35*1i*omega+1);
% % 
% % plot(20*log10(omega), 20*log10(abs(g2))); % Conversione in decibel
% % xlabel('Frequency (rad/s)');
% % ylabel('Magnitude (dB)');
% % title('Bode Plot of the 0D G2 magnitude');
% % grid on;
% 
% figure;
% plot(20*log10(omega), 20*log10(abs(g1))); % Conversione in decibel
% xlabel('Frequency (rad/s)');
% ylabel('Magnitude (dB)');
% title('Bode Plot of the 0D G1 magnitude');
% grid on;