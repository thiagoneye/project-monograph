%{
Federal University of Paraiba
Technology Center
Mechanical Engineering Department

Undergraduate Thesis
 
Numerical Analysis of Transient Conduction with Variable Source Term in
Fuel Rods of Nuclear Reactors by the Line Method 

Student: Thiago Ney Evaristo Rodrigues
Advisor: Dr. Jacques Cesar dos Santos
%}

clear
close all
clc

%% Inputs

nr = 400; % Points in spatial grid
nt = 400; % Points in temporal grid
r0 = 0.001; % Beginning of the r axis
rl = 1; % End of the r axis (Length L)
t0 = 0; % Start time
tl = 1.5; % End time
Bi = 15; % Biot Number
Gast = 32.4;

%% Calculations

r = linspace(r0,rl,nr); % r axis
t = linspace(t0,tl,nt); % Time
theta0 = Gast*(1 - r.^2)/4 + Gast/(2*Bi);     % Initial temperature

[~,theta] = ode15s(@d2udx2,t,theta0); % Solution

%% Plots

[R,Y] = meshgrid(t,r);

figure
surf(R,Y,theta', 'EdgeColor', 'none')
xlabel('Time [s]')
ylabel('Radius')
zlabel('Temperature')

figure
surf(R,Y,theta', 'EdgeColor', 'none')
xlabel('Time [s]')
ylabel('Radius')
zlabel('Temperature')
view(0,90)
hcb = colorbar;
hcb.Title.String = "T";

vt = [0; 0.5; 1; 1.5];
v = zeros(length(vt),1);

for i = 1:length(vt)
    
    v(i) = find(t == vt(i));
    
end

figure
plot(r,theta(v(1),:), r,theta(v(2),:), r,theta(v(3),:), r,theta(v(4),:))
legend('t = 0', 't = 0.5', 't = 1', 't = 1.5')
xlabel('Radius')
ylabel('Temperature')
grid