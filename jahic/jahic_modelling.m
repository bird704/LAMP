clear; close all; clc;

% model development(figure 3)

% important variables : 
% X = concentration of biomass (g/L)
% S = concentration of substrate (g/L)
% OCR = oxygen consumption rate (g/h)

%% description for the vector components
% p(1) = Cs         % y(1) = X          
% p(2) = Cx         % y(2) = S      
% p(3) = Ks         % y(3) = OCR
% p(4) = qm
% p(5) = qs_max
% p(6) = Si
% p(7) = Yem
% p(8) = Yos_an
% p(9) = Yos_en
                         
% set time step
    % t1 = glycerol batch
    % t2 = glycerol fed-batch
    % t3 = methanol fed-batch
    % t4 = constant methanol fed-batch

step = 10;
split = [25, 27.5, 30.7, 32.7, 40];
t1 = linspace(split(1),split(2),step);
t2 = linspace(split(2),split(3),step);
t3 = linspace(split(3),split(4),step);
t4 = linspace(split(4),split(5),step);

% Set parameter values (Jahic et al., 2002)
    % p(1) = Cs (carbon concentration in the substrate, g/g)     
    % p(2) = Cx (carbon cncentration in biomass, g/g)
    % p(3) = Ks (saturation constant, g/L)
    % p(4) = qm (maintainence coeff, g/(g*h)) 
    % p(5) = qs_max (specific max rate of substrate consumption, g/(g*h))
    % p(6) = Si (inlet substrate concentration, g/L)
    % p(7) = Yem (biomass yield coeff exclusive maintenance, g/g)
    % p(8) = Yos_an (coeff for o2 consump per substrate for ana, g/g)
    % p(9) = Yos_en (coeff for o2 consump per substrate for energy, g/g)

p1 = [0.391,0.396,0.1,0,0.37,555,0.7,0,1.217];
p2 = [0.391,0.396,0.1,0,0.37,555,0.7,0,1.217];
p3 = [0.375,0.96,0.1,0.013,0.57,780,0.36,0.5,1.5];
p4 = [0.375,0.96,0.1,0.013,0.57,780,0.36,0.5,1.5];

%% period 1(gly batch)
% Set initial values

% y_01(1) = X(Initial biomass conc., g/L)
% y_01(2) = S(Initial substrate conc. g/L)
% y_01(3) = OCR(Initial oxygen consumption rate g/h)

y_01 = [15,16,6.6];

% solve systems of ODE(for glycerol period)
[t1, y1] = ode15s(@(t,y) ODEsystem_jahic(t, y, p1, split), t1, y_01);

%% period 2(gly fed-batch)
y_02 = y1(step,:);

% solve systems of ODE(for methanol period)
[t2, y2] = ode15s(@(t,y) ODEsystem_jahic(t, y, p2, split), t2, y_02);

%% period 3(meth fed-batch)
y_03 = y2(step,:);

% solve systems of ODE(for methanol period)
[t3, y3] = ode15s(@(t,y) ODEsystem_jahic(t, y, p3, split), t3, y_03);

%% period 4(constant meth fed-batch)
y_04 = y3(step,:);

% solve systems of ODE(for methanol period)
[t4, y4] = ode15s(@(t,y) ODEsystem_jahic(t, y, p4, split), t4, y_04);

% plotting

t_all = vertcat(t1,t2,t3,t4);
y_all = vertcat(y1,y2,y3,y4);

subplot(2,3,1)
plot(t_all, y_all(:,1)); 
xlabel('t')
ylabel('g/L')
legend('X')

subplot(2,3,2)
plot(t_all, y_all(:,2)); 
xlabel('t')
ylabel('g/L')
legend('S')

subplot(2,3,3)
plot(t_all, y_all(:,3)); 
xlabel('t')
ylabel('g/h')
legend('OCR')

subplot(2,3,4)
plot(t_all, y_all(:,1:3)); 
xlabel('t')
legend('X','S','OCR')

%% DEBUGGING

t_check = linspace(25,40,100);
F_in = zeros(length(t_check),1);

for i = 1:length(t_check)
    F_in(i) = feed_rate(t_check(i), split)
end

subplot(2,3,5)
plot(t_check, F_in); 
xlabel('t')
ylabel('L/h')
legend('F_in')




