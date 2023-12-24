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
tl = 1.5;    % End time
Bi = 40;     % Biot Number
Gast = 32.4; % Heat source term

%% Curves

curve01 = [0.00609	9.10682;
0.09589	9.00000;
0.18722	8.83976;
0.30441	8.35905;
0.39422	7.87834;
0.50837	7.02374;
0.59209	6.27596;
0.70167	5.12760;
0.808219	3.81899;
0.910198	2.37685;
0.984779	1.25519];

curve02 = [0.00000	9.37389;
0.09741	9.26706;
0.21918	9.05341;
0.30898	8.73294;
0.39878	8.30564;
0.51750	7.58457;
0.60274	6.89021;
0.71385	5.79525;
0.802131	4.7003;
0.890411	3.25816;
0.998478	1.38872];

curve03 = [0.00000	11.91100;
0.09741	11.80420;
0.18265	11.67060;
0.33181	11.10980;
0.39878	10.73590;
0.52816	9.72107;
0.62557	8.70623;
0.69406	7.85163;
0.79300	6.30267;
0.90563	3.95252;
0.99696	1.81602];

curve04 = [0.00000	15.70330;
0.10198	15.56970;
0.20396	15.30270;
0.28767	14.87540;
0.39726	14.10090;
0.51294	12.95250;
0.59817	11.80420;
0.70929	9.98813;
0.799087	8.09199;
0.893455	5.50148;
0.996956	2.35015];

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
    Bi = 40;
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
    c2 = 1;
    c3 = 1;
    Glin = Gast*(1 + c2*r.^2)*exp(c3*t);
    
    for i = 1:nr
        ut(i) = (1/r(i))*ur(i) + urr(i) + Glin(i);
    end
end