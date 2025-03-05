%Registro de entrada
g_r=load("audiovoz05.txt", "-ascii");

%Parametros de entrada
N=length(g_r);
fs=16000;
dt=1/fs;
t0=0;
Tp=N/fs;
dw=(2*pi)/Tp;

disp("Dt:"); dt
disp("\nN:"); N
disp("\nDw:"); dw

%Vector tiempo
for i = 1:N
    t(i) = dt*(i-1);
 endfor
%Grafica
figure(1)
 plot(t, g_r, 'r')
 %TDF
 G_r_K=fft(g_r,N);

G_r_K_mod=abs(G_r_K(1:N/2));

figure(2);
stem(G_r_K_mod(1:2500), "b");
title("Módulo de G(k)");

%FILTRO
p=200*dw;

for i = 1:N
    h_1(i) = e^(-p*t(i));
  endfor

for i = 1:N
    h_2(i) = t(i)*(e^(-p*t(i)));
  endfor

h_t=dt*conv(h_1,h_2);

figure(3);
plot(t(1:320), h_t(1:320), 'r');
title('Respuesta al impulso unitario del filtro');

H_K=fft(h_t,N);

H_K_mod=abs(H_K(1:N/2));

figure(4);
stem(H_K_mod(1:500), "b");
title("Módulo de H(K)");

y_f=dt*conv(h_t, g_r);
Y_F=fft(y_f,N);
Y_F_mod=abs(Y_F(1:N/2));

figure(5)
plot(t(1:N), y_f(1:N));

figure(6);
stem(Y_F_mod(1:2500), "b");
title("Módulo de Y_F(k)");

