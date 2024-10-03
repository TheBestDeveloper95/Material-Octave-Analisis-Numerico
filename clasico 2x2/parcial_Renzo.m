function parcialote_Renzo

    %Definiciones
    y10=0;     % Valor Inicial VECTOR SOLUCION de funcion x1(t)
    y20=0;     % Valor Inicial VECTOR SOLUCION de funcion x2(t)
    t0 = 0.000;   %valor inicial tiempo
    tf = 12.5664;   %valor final en el tiempo
    w=1/2;

    N = 3000;   %cantidad de incrementos a realizar
    dt = (tf - t0)/N; %incremento en el tiempo
    dim = N + 1;

    y = zeros(2, N);     % Matriz para el vector solucion y(t)
    x = zeros(1, N);
    t=zeros(1,N);   % vector fila para el tiempo

    %inicializaciones
    t(1)=t0;
    y(1,1)=y10;
    y(2,1)=y20;


    %Euler EXPLÍCITO (Adelante)

    for i=1 : dim - 1      %((((Debería ser N-1?))))
        k1 = dt * f_pend(t(i),y(:,i));          %Euler IMPLÍCITO (Atrás)      k1 = dt * f_pend(t(i+1),y(:,i+1));
        y(:,i+1) = y(:,i) + k1;
        t(i+1) = t(i) + dt;
    endfor

    figure(1)
    plot(t, y(1,:), "black");
    grid on
    legend('x1(t)')

    figure(2)
    plot(t, y(2,:), "black");
    grid on
    legend('x2(t)')


#{

% Dimensionamiento RUNGE KUTTA
k1=zeros(2,1);
k2=zeros(2,1);

% RUNGE KUTTA

for j=1:N-1         %((((Debería ser N-1?))))

   k1 = dt * f_pend(t(j) , y(:,j));

   tg = t(j) + (dt/(2*w));
   yg = y(:,j) + (k1)/(2*w);

   k2 = dt * f_pend(tg,yg);

   y(:,j+1) = y(:,j) + (1-w) * k1 + w * k2;
   t(j+1)    = t(j) + dt;

end

 #}



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


figure(3)
plot(t,dd,'-r')
grid on
legend('Derivada de x1(t)')

    %Valores concretos pedidos sobre la funcion

valPed=2;
jota=N/valPed;

display("El valor de x1(t1) es :"),display(y(1, round(jota)));
display("El valor de x2(t1) es :"),display(y(2, round(jota)));
display("En t1=2seg el j es :"),display(round(jota));
display("En j=10001 t es :"),display(t(1,round(jota)));


    %Integral de Y1*Y1

    I11 = 0;

    for i=1 : dim
        func1(i) = y(1,i) * y(1,i);
        if(i >= 2)
            I11 = I11 + ((func1(i-1) + func1(i))/2) * dt;
        endif
    endfor

display("El valor de la integral definida () es :"),display(I11);

    %Integral I12

    I12 = 0;

    for i=1 : dim
        func2(i) = y(1,i) * dd(1,i);
        if(i >= 2)
            I12 = I12 + ((func2(i-1) + func2(i))/2) * dt;
        endif
    endfor

display("El valor de la integral definida () es :"),display(I12);

    % I11/I12 ver si es I12/I11

display("El valor del cociente es :"),display(I12/I11);

    %Método simpson

    ISimp = 0;

    for i=1:2:dim-1
        ISimp = ISimp + (y(1,i) + 4*y(1,i+1) + y(1,i+2)) * (dt/3);
    endfor

    display(ISimp);

    %Comparación con el de los trapecios para corroborar que esté bien

    ITrap = 0;

    for i=1:dim-1
        ITrap = ITrap + ((y(1,i) + y(1,i+1))/2) * dt;
    endfor

    display(ITrap);

end

function [fy] = f_pend(x,z)

  fy(1,1)=(0*z(1) + 1*z(2) )+ 0*sin(3*x);
  fy(2,1)=(-((2)^2)*z(1) - 0.020*z(2)) - 1*sin(3*x);

end
