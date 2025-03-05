function intento2
  clear;
  clc;

  % Cargando el registro de datos
  g_r = load("audiovoz02.txt", "-ascii");

  N = length(g_r); fs = 16000; dt = 1 / fs;
  t0 = 0; Tp = N / fs; dw = (2 * pi) / Tp;

  disp("Dt:"); dt
  disp("\nN:"); N
  disp("\nDw:"); dw

  %%===========================ARREGLOS LUEGO UTILIZADOS========================%%

  t = (0:N-1) * dt;

  % Grafica del registro
  figure(1);
  plot(t, g_r, 'r');
  title("Señal Original");

  %=================================TDF=========================================%%
  G_r_K = fft(g_r);

  G_r_K_mod = abs(G_r_K(1:N/2));

  figure(2);
  stem(G_r_K_mod, "b");
  title("Módulo de G(k)");

  Adm = max(G_r_K_mod);
  i_max = find(G_r_K_mod == Adm); % índice correspondiente a la máxima amplitud.
  fm = (i_max - 1) * dw; % frecuencia angular correspondiente. El menos 1 es por que las frecuencias arrancan desde 0 pero los arreglos desde 1

  disp("\nAdm:"); printf('%.2f\n', Adm);
  disp("\nfm:"); printf('%.2f\n', fm);

  %==============================FILTROS======================================%%
  zita = 0.4;
  wn = fm;

  % FUNCIÓN ESCALÓN
  u = zeros(N, 1);
  t0 = Tp / 2; % Tiempo en el que el escalón se activa
  for i = 1:N
      if t(i) < t0
          u(i) = 0; % Antes del umbral
      else
          u(i) = 1; % Después del umbral
      end
  end

  % Resolución del sistema con Euler explícito
  x3_2 = zeros(3, N);
  for i = 1:N-1
      k1 = dt * f_pend(u(i), x3_2(:,i), wn, zita); % Paso de Euler con la función escalón como entrada
      x3_2(:,i+1) = x3_2(:,i) + k1;
  end

  y_u = x3_2(3, :); % Salida del sistema

  %===============================OBTENCION DE H================================%%
  h = zeros(N, 1); % Inicializar el vector de la derivada

  % DERIVADA PRIMERA ASIMÉTRICA (inicio)
  h(1) = (1 / (2 * dt)) * (-3 * y_u(1) + 4 * y_u(2) - y_u(3));

  % DERIVADA PRIMERA CENTRAL (para puntos intermedios)
  for i = 2:(N-1)
      h(i) = (1 / (2 * dt)) * (y_u(i+1) - y_u(i-1));
  end

  % TDF de h(t)
  H_k = fft(h);

  H_k_mod = abs(H_k(1:N/2));

  Ah = H_k_mod(i_max);

  disp("\nAh:"); Ah

  % Graficando

   figure(3)
  subplot(2,1,1), plot(t,h, "b" )
  grid on
  title ('Función h(t)')
  subplot(2,1,2), stem(H_k_mod, "b")
  grid on
  title ('Módulo de H(k)')

 %Reinicio de los arreglos
x_g = zeros(3, N); % Estados iniciales del sistema

for i = 1:N-1
    k1 = dt * f_pend(g_r(i), x_g(:,i), wn, zita); % Paso de Euler con u = g
    x_g(:,i+1) = x_g(:,i) + k1;
end

% Resultado del sistema con g como entrada
y_f = x_g(3, :); % Salida del sistema para u = g

  % Transformada discreta de Fourier de la salida
  Y_F = fft(y_f);
  Y_F_mod = abs(Y_F(1:N/2));


   figure(4)
  subplot(2,1,1), plot(t,y_f(1:N), "b" )
  grid on
  title ('Salida del sistema con convolución')
  subplot(2,1,2), stem(Y_F_mod, "b")
  grid on
  title ('Módulo de Y(k)')

  Ay = Y_F_mod(i_max);
  disp("\nAy: "); Ay
endfunction

function [fy] = f_pend(u, z, wn, zita)
  fy = zeros(3, 1);
  fy(1) = z(2) + u;
  fy(2) = -wn^2 * z(1) - 2 * zita * wn * z(2);
  fy(3) = z(1) - wn * z(3);
end
