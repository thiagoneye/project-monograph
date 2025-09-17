% UNIVERSIDADE FEDERAL DA PARAÍBA
% CENTRO DE TECNOLOGIA
% DEPARTAMENTO DE ENGENHARIA MECÂNICA
%
% DISCENTE THIAGO NEY EVARISTO RODRIGUES
% ORIENTADOR DR. JACQUES CÉSAR DOS SANTOS
%
% TRABALHO DE CONCLUSÃO DE CURSO
% 
% ANÁLISE NUMÉRICA DE CONDUÇÃO TRANSIENTE COM TERMO FONTE VARIÁVEL EM
% VARETAS COMBUSTÍVEIS DE REATORES NUCLEARES PELO MÉTODO DAS LINHAS 

clear
close all
clc

%% Inputs

nr = 100;    % Points in spatial grid
nt = 100;    % Points in temporal grid
r0 = 0.001;  % Beginning of the r axis
rl = 1;      % End of the r axis (Length L)
t0 = 0;      % Start time
tl = 1.5;    % End time
Bi = 15;     % Biot number
Gast = 32.4; % Heat source term

%% Curves

curve01 = [0.00174	9.22535;
0.11130	9.08451;
0.19478	8.92019;
0.29739	8.49765;
0.40174	7.91080;
0.50957	7.11268;
0.58957	6.38498;
0.709565 5.14085;
0.808696 3.92019;
0.909565 2.48826;
0.998261 1.17371];

curve02 = [0.00348 3.61502;
0.08174	3.59155;
0.20348	3.52113;
0.28348	3.40376;
0.40348	3.14554;
0.48348	2.93427;
0.60348	2.53521;
0.68174	2.23005;
0.80000	1.69014;
0.91652	1.00939;
0.99130	0.53991];

curve03 = [0.00174	7.30047;
0.08174	7.25352;
0.20174	7.06573;
0.28174	6.83099;
0.40000	6.36150;
0.51652	5.70423;
0.59130	5.18779;
0.70087	4.29577;
0.80696	3.26291;
0.90783	2.08920;
0.972174 1.24413];

curve04 = [0.00348	11.80750;
0.12174	11.64320;
0.20174	11.43190;
0.31826	10.82160;
0.39304	10.28170;
0.49913	9.34272;
0.60348	8.19249;
0.69739	6.90141;
0.78783	5.49296;
0.89913	3.47418;
0.97913	1.90141];


%% Calculations

r = linspace(r0,rl,nr);                   % r axis
t = linspace(t0,tl,nt);                   % Time
theta0 = Gast*(1 - r.^2)/4 + Gast/(2*Bi); % Initial temperature
[~,theta] = ode15s(@d2udx2,t,theta0);     % Solution

%% Plots

[R,Y] = meshgrid(t,r);

figure
surf(R,Y,theta', 'EdgeColor', 'none')
xlabel('Tempo [s]')
ylabel('Raio')
zlabel('Temperatura')

figure
surf(R,Y,theta', 'EdgeColor', 'none')
xlabel('Tempo [s]')
ylabel('Raio')
zlabel('Temperatura')
view(0,90)
hcb = colorbar;
hcb.Title.String = "T";

vt = [0; 0.5; 1; 1.5];
v = zeros(length(vt),1);

for i = 1:length(vt)
    
    v(i) = find(t == vt(i));
    
end

figure
plot(r,theta(v(1),:), 'b', 'Color', "#EDB120")
hold on
plot(curve01(:,1), curve01(:,2), 'x', 'Color', "#EDB120")
plot(r,theta(v(2),:), 'g', 'Color', "#0072BD")
plot(curve02(:,1), curve02(:,2), 'xg', 'Color', "#0072BD")
plot(r,theta(v(3),:), 'r', 'Color', "#D95319")
plot(curve03(:,1), curve03(:,2), 'xr', 'Color', "#D95319")
plot(r,theta(v(4),:), 'y', 'Color', "#7E2F8E")
plot(curve04(:,1), curve04(:,2), 'xy', 'Color', "#7E2F8E")
hold off
legend('t = 0', 't = 0 (Soares2017)', ...
    't = 0.5', 't = 0.5 (Soares2017)', ...
    't = 1', 't = 1 (Soares2017)', ...
    't = 1.5', 't = 1.5 (Soares2017)')
xlabel('Raio')
ylabel('Temperatura')
grid

%% Function

function ut = d2udx2(t,u)
	% Problem parameters
	r0 = 0.001;
    rl = 1;
    Bi = 15;
    Gast = 32.4;
    nr = length(u);
    r = linspace(r0,rl,nr);
    
    % Partial derivative
	dr = (rl - r0)/(nr-1);
    dr2 = dr^2;
    
    ur = zeros(nr,1); % Preallocation
	urr = zeros(nr,1); % Preallocation
	for i = 2:(nr-1)
        ur(i) = (u(i) - u(i-1))/dr;
        urr(i) = (u(i+1) - 2*u(i) + u(i-1))/dr2;
	end

	% Boundary condition
    ur(1) = 0;
	urr(1) = 2*(u(2) - u(1))/dr2;
    ur(end) = (u(end) - u(end - 1))/dr;
	urr(end) = 2*(u(end - 1) - (Bi*dr + 1)*u(end))/dr2;

    % Partial differential equation
    ut = zeros(nr,1);
    Glin = Gast*t;
    
    ut(1) = 2*urr(1) + Glin;
    for i = 2:nr
        ut(i) = (1/r(i))*ur(i) + urr(i) + Glin;
    end
  
end