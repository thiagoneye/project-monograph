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

r0 = 0.001;  % Beginning of the r axis
rl = 1;      % End of the r axis (Length L)
t0 = 0;      % Start time
tl = 1.5;    % End time
Bi = 15;     % Biot number
Gast = 32.4; % Heat source term

nr = 86;                           % Points in spatial mesh
number_of_meshes = 13;             % Number of temporal meshes
mesh = zeros(number_of_meshes, 1); % Preallocation of the vector of the number of points in the temporal mesh

%% Calculations

% Initial Parameters
mesh(1) = 11;                                              % Initial number of points in the temporal mesh
convergence_matrix = zeros(mesh(1), nr, number_of_meshes); % Convergence matrix

for idx = 2:number_of_meshes
    mesh(idx) = (mesh(idx-1) - 1)*2 + 1;
end

% Simulation
for idx = 1:number_of_meshes
    nt = mesh(idx);                           % Points in temporal mesh
    r = linspace(r0,rl,nr);                   % r axis
    t = linspace(t0,tl,nt);                   % Time
    theta0 = Gast*(1 - r.^2)/4 + Gast/(2*Bi); % Initial temperature
    
    [~,theta] = ode15s(@d2udx2,t,theta0);     % Solution

    index_vector = linspace(1, mesh(idx), mesh(1));
    convergence_matrix(:, :, idx) = theta(index_vector, :);
end

% Error
error_matrix = diff(convergence_matrix, 1, 3);
maximum_error = reshape(max(error_matrix, [], [1 2]), size(error_matrix, 3), 1);

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
    Glin = Gast*t;
    for i = 1:nr
        ut(i) = (1/r(i))*ur(i) + urr(i) + Glin;
    end
  
end