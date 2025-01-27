% Parámetros 
prmtr.Ip = 0.0079;      % kg*m^2  Inerciap
prmtr.Mc = 0.7031;      % kg      Masac
prmtr.lp = 0.3302;      % m       Dist del pivote al CG del p
prmtr.Mp = 0.23;        % kg      Masap
prmtr.Fc = 0;           % N       Fuerza externa 
prmtr.Beq = 4.3;        % Ns/m    Coef de amortiguamiento del c
prmtr.g = 9.81;         % m/s^2   Aceleración gravitacional
prmtr.Bp = 0.0024;      % Nms/rad Coef de amortiguamiento del p

% Condiciones iniciales
alpha0 = deg2rad(1);  % grados a radianes
xc0 = 0;              % Posición inicial c
alphadot0 = 0;        % Velocidad angular inicial p
xcdot0 = 0;           % Velocidad inicial c

x0 = [xc0; xcdot0; alpha0; alphadot0]; % Vector de condiciones iniciales

% Tiempo de simulación
tspan = [0 20];

% ode45 
[t, solu] = ode45(@(t, x) [
    x(2); 
    ((prmtr.Ip + prmtr.Mp*prmtr.lp^2)*prmtr.Fc + prmtr.Mp^2*prmtr.lp^2*prmtr.g*cos(x(3))*sin(x(3)) ...
    - (prmtr.Ip + prmtr.Mp*prmtr.lp^2)*prmtr.Beq*x(2) - (prmtr.Ip*prmtr.Mp*prmtr.lp - prmtr.Mp^2*prmtr.lp^3)*x(4)^2*sin(x(3)) ...
    - prmtr.Mp*prmtr.lp*x(4)*cos(x(3))*prmtr.Bp) / ...
    ((prmtr.Mc + prmtr.Mp)*prmtr.Ip + prmtr.Mc*prmtr.Mp*prmtr.lp^2 + prmtr.Mp^2*prmtr.lp^2*sin(x(3))^2);
    x(4);
    ((prmtr.Mc + prmtr.Mp)*prmtr.Mp*prmtr.g*prmtr.lp*sin(x(3)) - (prmtr.Mc + prmtr.Mp)*prmtr.Bp*x(4) ...
    + prmtr.Fc*prmtr.Mp*prmtr.lp*cos(x(3)) - prmtr.Mp^2*prmtr.lp^2*x(4)^2*sin(x(3))*cos(x(3)) ...
    - prmtr.Beq*prmtr.Mp*prmtr.lp*x(2)*cos(x(3))) / ...
    ((prmtr.Mc + prmtr.Mp)*prmtr.Ip + prmtr.Mc*prmtr.Mp*prmtr.lp^2 + prmtr.Mp^2*prmtr.lp^2*sin(x(3))^2)
], tspan, x0);

% Graficar resultados

% Desplazamiento del carrito
figure;
subplot(2,2,1);
plot(t, solu(:,1), 'Color', [1 0.4 0.6], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Desplazamiento xc (m)');
title('Desplazamiento del carrito');
grid on;

% Ángulo del péndulo
subplot(2,2,2);
plot(t, solu(:,3), 'Color', [1 0.6 0.8], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Ángulo α (rad)');
title('Ángulo del péndulo');
grid on;

% Velocidad del carrito (xc_dot)
subplot(2,2,3);
plot(t, solu(:,2), 'Color', [1 0.2 0.6], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Velocidad xc˙ (m/s)');
title('Velocidad del carrito');
grid on;

% Velocidad angular del péndulo (alpha_dot)
subplot(2,2,4);
plot(t, solu(:,4), 'Color', [0.8 0.2 0.5], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Velocidad angular α˙ (rad/s)');
title('Velocidad angular del péndulo');
grid on;
