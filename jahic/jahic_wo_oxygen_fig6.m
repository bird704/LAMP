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
muu = 0.55;
split = [0, 27.5, 31, 32.5, 34, 100];
t1 = linspace(split(1),split(2),step);
t2 = linspace(split(2),split(3),step);
t3 = linspace(split(3),split(4),step);
t4 = linspace(split(4),split(5),step);
t5 = linspace(split(5),split(6),step);

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

qm_vec = [0, 0.013, 0.04];
for i=1:3
    p1 = [0.391,0.396,0.1,0,0.37,555,0.7,0,1.217];
    p2 = [0.391,0.396,0.1,0,0.37,555,0.7,0,1.217];
    p3 = [0.375,0.96,0.1,qm_vec(i),0.57,780,0.36,0.5,1.5];
    p4 = [0.375,0.96,0.1,qm_vec(i),0.57,780,0.36,0.5,1.5];
    p5 = [0.375,0.96,0.1,qm_vec(i),0.57,780,0.36,0.5,1.5];
%% period 1(gly batch)
% Set initial values
% y_01(1) = X(Initial biomass conc., g/L)
% y_01(2) = S(Initial substrate conc. g/L)
% y_01(3) = OCR(Initial oxygen consumption rate g/h)
    y_01 = [0.025,38,0];  % X0 = 0이면 시뮬레이션이 안돌아가는 문제?

% solve systems of ODE(for glycerol period)
    [t1, y1] = ode15s(@(t,y) ODEsystem_jahic(t, y, p1, split,muu), t1, y_01);

%% period 2(gly fed-batch)
    y_02 = y1(step,:);

% solve systems of ODE(for methanol period)
    [t2, y2] = ode15s(@(t,y) ODEsystem_jahic(t, y, p2, split,muu), t2, y_02);

%% period 3(meth fed-batch 1)
    y_03 = y2(step,:);

% solve systems of ODE(for methanol period)
    [t3, y3] = ode15s(@(t,y) ODEsystem_jahic(t, y, p3, split,muu), t3, y_03);

%% period 4(meth fed-batch 2)
    y_04 = y3(step,:);

% solve systems of ODE(for methanol period)
    [t4, y4] = ode15s(@(t,y) ODEsystem_jahic(t, y, p4, split,muu), t4, y_04);

%% period 5(constant meth fed-batch)
    y_05 = y4(step,:);

% solve systems of ODE(for methanol period)
    [t5, y5] = ode15s(@(t,y) ODEsystem_jahic(t, y, p5, split,muu), t5, y_05);

    t_all = vertcat(t1,t2,t3,t4,t5);
    y_all = vertcat(y1,y2,y3,y4,y5);
    
    hold on
    subplot(2,2,1)
    plot(t_all, y_all(:,1)); 
    xlabel('t'); ylabel('g/L')
    ylim([0,200])

    hold on
    subplot(2,2,2)
    plot(t_all, y_all(:,2)); 
    xlabel('t'); ylabel('g/L')
    ylim([0,40])

    hold on
    subplot(2,2,3)
    plot(t_all, y_all(:,3)); 
    xlabel('t'); ylabel('g/h')
    ylim([0,32])

end

hold off

subplot(2,2,1)
legend('X, qm=0', 'X, qm=0.013', 'X, qm=0.04')
subplot(2,2,2)
legend('S')
subplot(2,2,3)
legend('OCR, qm=0', 'OCR, qm=0.013', 'OCR, qm=0.04')


%% DEBUGGING F_in

t_check = linspace(0,100,500);
F_in = zeros(length(t_check),1);

for i = 1:length(t_check)
    F_in(i) = feed_rate(t_check(i), split, muu);
end

subplot(2,2,4)
plot(t_check, F_in); 
xlabel('t')
ylabel('L/h')
legend('F_in')





