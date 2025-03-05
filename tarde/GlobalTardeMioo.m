function Globalazo

format bank;
 %-----------      Generación de la función discreta dato     -----------
  gr=load("audiovoz05.txt", "-ascii");

N=length(gr); fs=16000; dt=1/fs;
ti=0; Tp=N/fs; dw=(2*pi)/Tp;

 for i = 1:N
    tn(i) = ti + dt*(i-1);
  endfor

disp("\nEl delta tiempo 'Dt' es:"); printf('%.8f\n', dt);
disp("\nLa cantidad de elementos 'N' es:"); N
disp("\nEl delta de frecuencia angular 'Δw' es:"); printf('%.3f\n', dw);


    %-----------      Transformada discreta de Fourier      -----------

gr_tdf=fft(gr);
g_tdf_mod=abs(gr_tdf(1:N));

  figure(1)
  subplot(2,1,1), plot(tn,gr, "b" )
  grid on
  title ('gr(t) Entrada')
  subplot(2,1,2), stem(g_tdf_mod(1:2500), "b")
  grid on
  title ('Módulo de TDF de gr(t) Entrada')

                %-----------      Convolución Parte 1    -----------

p= 200 * dw;
  for i = 1:N
  h1(i)= e^(-p*tn(i));
endfor
  for i = 1:N
  h2(i)= tn(i) * e^(-p*tn(i));
endfor

h= conv(h1,h2)*dt;
h_tdf=fft(h,N);
h_tdf_mod=abs(h_tdf(1:N/2));

figure(2)
  subplot(2,1,1), plot( tn(1:320), h(1:320), "r" ) %fs/0.02
  grid on
  title ('h(t) Respuesta a Impulso Unitario del Filtro')
  subplot(2,1,2), stem(h_tdf_mod(1:500), "b")
  grid on
  title ('Módulo de TDF de h(t)')

                %-----------      Convolución Parte 2    -----------

  yf_= dt*conv(h , gr);
  Y_tdf= fft(yf_,N);
  yf_tdf_mod=abs(Y_tdf(1:N/2));

 figure(3)
  subplot(2,1,1), plot( tn(1:N),yf_(1:N), "r" )
  grid on
  title ('h(t) Respuesta a Impulso Unitario del Filtro')
  subplot(2,1,2), stem(yf_tdf_mod(1:2500), "r")
  grid on
  title ('Módulo de TDF de h(t)')
endfunction
