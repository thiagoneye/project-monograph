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
nt = 301;    % Points in temporal grid
r0 = 0;      % Beginning of the r axis
rl = 1;      % End of the r axis (Length L)
t0 = 0;      % Start time
tl = 0.5;    % End time
Bi = 15;     % Biot Number
Gast = 32.4; % Heat source term

%% Curves

curve01 = [0.00609	9.22656;
0.09589	9.11830;
0.21461	8.84952;
0.30289	8.47462;
0.41400	7.81486;
0.49620	7.20896;
0.59970	6.28267;
0.69863	5.2498;
0.791476	4.11037;
0.896499	2.72184;
0.987823	1.29799];

curve02 = [0.00000	7.75111;
0.11263	7.64244;
0.18722	7.51667;
0.29985	7.17689;
0.40792	6.69497;
0.51294	6.07088;
0.61339	5.32243;
0.71233	4.46734;
0.80670	3.50566;
0.891933	2.43748;
0.977169	1.31596];

curve03 = [0.00000	5.51111;
0.11263	5.50911;
0.18722	5.50778;
0.30137	5.48798;
0.41400	5.41486;
0.49011	5.32462;
0.60122	5.03820;
0.71081	4.55625;
0.81126	3.80780;
0.908021	2.73959;
0.980213	1.5648];

curve04 = [0.00152	9.31553;
0.11263	9.43800;
0.18722	9.56112;
0.29833	9.89692;
0.40944	10.23270;
0.48402	10.42700;
0.59665	10.42490;
0.70472	9.92525;
0.79148	8.89260;
0.91411	5.95726;
0.946073	4.94336];

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

vt = [0; 0.05; 0.25; 0.5];
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
    't = 0.05', 't = 0.05 (Soares2017)', ...
    't = 0.24', 't = 0.24 (Soares2017)', ...
    't = 0.5', 't = 0.5 (Soares2017)')
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
    ur(1) = (u(1) - u(2))/dr;
	urr(1) = 2*(u(2) - u(1))/dr2;
    ur(end) = (u(end) - u(end - 1))/dr;
	urr(end) = 2*(u(end - 1) - (Bi*dr + 1)*u(end))/dr2;

    % Partial differential equation
    ut = zeros(nr,1);
    c3 = 4;
    Glin = Gast*(r.^2)*exp(c3*t);
    
    for i = 1:nr
        ut(i) = (1/r(i))*ur(i) + urr(i) + Glin(i);
    end
end