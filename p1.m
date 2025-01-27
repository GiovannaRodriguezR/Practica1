% Parámetros del sistema
params.Ip = 0.0079;      % kg*m^2  Momento de inercia del péndulo
params.Mc = 0.7031;      % kg      Masa del carrito
params.lp = 0.3302;      % m       Distancia pivote a CG
params.Mp = 0.23;        % kg      Masa del péndulo
params.Fc = 0;           % N       Fuerza aplicada
params.Beq = 4.3;        % Ns/m    Coeficiente de amortiguamiento equivalente
params.g = 9.81;         % m/s^2   Aceleración gravitacional
params.Bp = 0.0024;      % Nms/rad Coeficiente de amortiguamiento en el eje del péndulo

% Condiciones iniciales
alpha0 = deg2rad(1);  % Convertimos de grados a radianes
xc0 = 0;
alphadot0 = 0;
xcdot0 = 0;

x0 = [xc0; xcdot0; alpha0; alphadot0]; % Vector de condiciones iniciales

% Tiempo de simulación
tspan = [0 20];

% Resolver el sistema con ode45, pasando los parámetros
[t, sol] = ode45(@(t, x) pendulumODE(t, x, params), tspan, x0);

% Graficar los resultados

% Desplazamiento del carrito
figure;
subplot(2,2,1);
plot(t, sol(:,1), 'Color', [1 0.4 0.6], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Desplazamiento xc (m)');
title('Desplazamiento del carrito');
grid on;

% Ángulo del péndulo
subplot(2,2,2);
plot(t, sol(:,3), 'Color', [1 0.6 0.8], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Ángulo α (rad)');
title('Ángulo del péndulo');
grid on;

% Velocidad del carrito (xc_dot)
subplot(2,2,3);
plot(t, sol(:,2), 'Color', [1 0.2 0.6], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Velocidad xc˙ (m/s)');
title('Velocidad del carrito');
grid on;

% Velocidad angular del péndulo (alpha_dot)
subplot(2,2,4);
plot(t, sol(:,4), 'Color', [0.8 0.2 0.5], 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Velocidad angular α˙ (rad/s)');
title('Velocidad angular del péndulo');
grid on;

% Función del sistema de ecuaciones diferenciales
function dx = pendulumODE(t, x, params)
    xc = x(1);
    xcdot = x(2);
    alpha = x(3);
    alphadot = x(4);

    % Extraer parámetros
    Ip = params.Ip;
    Mc = params.Mc;
    lp = params.lp;
    Mp = params.Mp;
    Fc = params.Fc;
    Beq = params.Beq;
    g = params.g;
    Bp = params.Bp;

    % Ecuaciones del sistema
    denominator = (Mc + Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(alpha)^2;
    
    xcdotdot = ((Ip + Mp*lp^2)*Fc + Mp^2*lp^2*g*cos(alpha)*sin(alpha) ...
               - (Ip + Mp*lp^2)*Beq*xcdot - (Ip*Mp*lp - Mp^2*lp^3)*alphadot^2*sin(alpha) ...
               - Mp*lp*alphadot*cos(alpha)*Bp) / denominator;
           
    alphadotdot = ((Mc + Mp)*Mp*g*lp*sin(alpha) - (Mc + Mp)*Bp*alphadot + ...
                   Fc*Mp*lp*cos(alpha) - Mp^2*lp^2*alphadot^2*sin(alpha)*cos(alpha) ...
                   - Beq*Mp*lp*xcdot*cos(alpha)) / denominator;
    
    dx = [xcdot; xcdotdot; alphadot; alphadotdot];
end
