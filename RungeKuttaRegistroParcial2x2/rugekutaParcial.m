function runge_kutta_sistemaRegistro_01

    %Definiciones
    y0_1 = 0;   % Valor Inicial VECTOR SOLUCION de funcion x1(t)
    y0_2 = 0;   % Valor Inicial VECTOR SOLUCION de funcion x2(t)
    t0 =   0;   % Valor inicial en el tiempo
    fs=5000;    % Frecuencia de muestreo
    dt = 0.0002;    %incremento en el tiempo

arreglo1d =load("p1_3k9_01.txt", "-ascii");
ndt=length(arreglo1d); %Defino N == longitud del arreglo

    dim= ndt;
    w = 0.5;    %Utilizo Runge Kutta (1/2) por la digitalización de la función g(t)

    t = zeros(1, ndt);    % vector para el tiempo
    y = zeros(2, ndt);    % Matriz para los vectores solucion y(t)
    x = zeros(2,ndt);
    k1 = zeros(2,1);
    k2 = zeros(2,1);

    %inicializaciones
    t(1)=t0;
    y(1,1) = y0_1;
    y(2,1) = y0_2;

    %RUNGE KUTTA

    for j=1: ndt-1  %N-1

        x=arreglo1d(j);
        k1 = dt * funcion(y(:,j), x);

        tg = t(j) + (dt/(2*w));
        yg = y(:,j) + k1/(2*w);

        x=arreglo1d(j+1);
        k2 = dt * funcion(yg, x);

        y(:, j+1) = y(:,j) + (1-w) * k1 + w * k2;
        t(j+1) = t(j) + dt;

    endfor

    figure(1)
    plot(t, y(1,:) , 'r')
    grid on
    legend('Función x1(t)')


        %Derivada central de orden 2 para x1(t)
%Declaraciones

for i=1:dim
  if (i==1)
  dd(i)=(-3*y(1,i)+4*y(1,i+1)-y(1,i+2))/(2*dt); %primera posicion
  elseif(i<dim)
  dd(i)=(-y(1,i-1)+y(1,i+1))/(2*dt); %desde i=2 hasta dim-1
  else
  dd(i)=(3*y(1,dim)-4*y(1,dim-1)+y(1,dim-2))/(2*dt);%ultima posicion
  endif
endfor


figure(2)
plot(t,dd,'-r')
grid on
legend('Derivada de x1(t)')

    %Valores concretos pedidos sobre la funcion

valPed=2;
jota=ndt/valPed;

display("El valor de x1(t1) es :"),display(y(1, round(jota)));
display("El valor de x2(t1) es :"),display(y(2, round(jota)));
display("En t1=2seg el j es :"),display(round(jota));
display("En j=10001 t es :"),display(t(1,round(jota)));
    %Integral

    Igg = 0;

    for i=1 : dim
        func1(i) = x * x;
        if(i >= 2)
            Igg = Igg + ((func1(i-1) + func1(i))/2) * dt;
        endif
    endfor

display("El valor de la integral definida (Igg) es :"),display(Igg);

    %Integral I11

    I11 = 0;

    for i=1 : dim
        func1(i) = y(1,i) * y(1,i);
        if(i >= 2)
            I11 = I11 + ((func1(i-1) + func1(i))/2) * dt;
        endif
    endfor

display("El valor de la integral definida (I11) es :"),display(I11);

    %Integral I22

    I22 = 0;

    for i=1 : dim
        func2(i) = dd(1,i) * dd(1,i);
        if(i >= 2)
            I22 = I22 + ((func2(i-1) + func2(i))/2) * dt;
        endif
    endfor

display("El valor de la integral definida (I22) es :"),display(I22);

      %Reproduzco
player=audioplayer(arreglo1d, fs, 8);
play(player);                              %Reproduzco audio
pause(5);

end


function [fy]=funcion(z,x)
    fy(1,1) = (0*z(1) + 1*z(2));
    fy(2,1) = (-14400*z(1) - 1.2*z(2))+ x;

end
