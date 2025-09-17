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

nr = 400;    % Points in spatial grid
nt = 400;    % Points in temporal grid
r0 = 0.001;  % Beginning of the r axis
rl = 1;      % End of the r axis (Length L)
t0 = 0;      % Start time
tl = 1.5;    % End time
Bi = 40;     % Biot Number
Gast = 32.4; % Heat source term

%% Curves

curve01 = [0.00159	9.19405;
0.12440	9.03756;
0.19777	8.84194;
0.29665	8.45070;
0.39394	7.90297;
0.50718	7.08138;
0.61563	6.06416;
0.71451	4.96870;
0.79904	3.91236;
0.905901	2.38654;
0.99681	0.978091];

curve02 = [0.00159	11.97180;
0.09729	11.81530;
0.20096	11.58060;
0.30463	11.07200;
0.39713	10.40690;
0.51675	9.23318;
0.60128	8.25509;
0.709729	6.76839;
0.805423	5.24257;
0.910686	3.16901;
0.995215	1.52582];

curve03 = [0.00000	16.39280;
0.09091	16.23630;
0.19936	15.88420;
0.31579	15.02350;
0.40829	14.04540;
0.50080	12.87170;
0.59330	11.42410;
0.69697	9.46792;
0.807018	7.04225;
0.909091	4.34272;
0.995215	1.99531;];

curve04 = [0.00159	20.97030;
0.10367	20.69640;
0.21053	20.18780;
0.29346	19.44440;
0.40510	17.99690;
0.50558	16.31460;
0.59490	14.51490;
0.701754	11.9327;
0.805423	8.99844;
0.896332	6.06416;
0.995215	2.58216;];

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
    ur(1) = 0;
	urr(1) = 2*(u(2) - u(1))/dr2;
    ur(end) = (u(end) - u(end - 1))/dr;
	urr(end) = 2*(u(end - 1) - (Bi*dr + 1)*u(end))/dr2;

    % Partial differential equation
    ut = zeros(nr,1);
    c1 = 1;
    Glin = Gast*(1 + c1*t);
    
    ut(1) = 2*urr(1) + Glin;
    for i = 2:nr
        ut(i) = (1/r(i))*ur(i) + urr(i) + Glin;
    end
  
end