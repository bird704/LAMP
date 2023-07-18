clear; close all; clc;

% model development

% important variables : 
% X = concentration of biomass
% S = concentration of substrate
% P = concentration of protein
% Cx = concentration of intracellular(X) carbon
% Nx = concentration of intracellular(X) nitrogen

%% description for the vector components
% p(1) = qm         % y(1) = X      % y(6)= ammonium(NT)      
% p(2) = Yem        % y(2) = S      % y(7)= glutamine(QT)
% p(3) = qs_max     % y(3) = P      % y(8)= arginine(RT)
% p(4) = Ks         % y(4) = Cx     % y(9)= phosphate(PT)
% p(5) = KaE        % y(5) = Nx     % y(10)= sulfate(ST)
% p(6) = H                          % y(11)= potassium(+)
% p(7) = a_P                        % y(12)= magnesium(2+)
% p(8) = a_C                        % y(13)= calcium(2+)
% p(9) = a_N                        % y(14)= trans M(2+)
                                    % y(15)= chloride(-)
% set time step

step = 10;
t1 = linspace(0,20,step);
t2 = linspace(20,40,step);

% Set parameter values

p1 = [0,0.7,0.37,0.1,10^(-3.42),10^(-6.5),0.0232,0,0];
p2 = [0.013,0.36,0.57,0.1,10^(-3.42),10^(-6.5),0.069,0.594,0.133];
    % p(1) = qm (Jahic et al., 2002)
    % p(2) = Yem (Jahic et al., 2002)
    % p(3) = qs_max (Jahic et al., 2002)
    % p(4) = Ks (Jahic et al., 2002)
    % p(5) = KaE (12. Hong moo sun)
    % p(6) = H   (exp condition)    
    % p(7) = a_P  (12. Hong moo sun)
    % p(8) = a_C
    % p(9) = a_N

    
% set initial values
y_01_p1 = [0 0 0 0.413 0.09];
y_01_p2 = [24.9 11.9 8.4 88.2 19.1 12.5 88.2 60.1 19.1 24.5 2.45];
y_01 = vertcat(transpose(y_01_p1),transpose(y_01_p2));

% y(6)= ammonium(NT)      
% y(7)= glutamine(QT)
% y(8)= arginine(RT)
% y(9)= phosphate(PT)
% y(10)= sulfate(ST)
% y(11)= potassium(+)
% y(12)= magnesium(2+)
% y(13)= calcium(2+)
% y(14)= trans M(2+)
% y(15)= chloride(-)


k1 = [3, 0.01, 555, 0.5, 0.2]; 
    % k1(1) = volume : 3 L   (Jahic et al., 2002) 
    % k1(2) = biomass(X)_in : zero?    
    % k1(3) = substrate(S)_in : 555 g/L glycerol  (matthews, 2017)
    % k1(4) = concentration of carbon in the protein
    % k1(5) = concentration of nitrogen in the protein


% solve systems of ODE(for glycerol period)
[t1, y1] = ode45(@(t,y) ODEsystem2(t, y, p1, k1), t1, y_01);



% solve systems of ODE(for methanol period)
y_02 = y1(step,:);
k2 = [3, 0.01, 780, 0.5, 0.2];
    % k2(1) = volume : 3 L   (Jahic et al., 2002) 
    % k2(2) = biomass(X)_in : zero?    
    % k2(3) = substrate(S)_in : 780 g/L glycerol  (matthews, 2017)
    % k2(4) = concentration of carbon in the protein
    % k2(5) = concentration of nitrogen in the protein

[t2, y2] = ode45(@(t,y) ODEsystem2(t, y, p2, k2), t2, y_02);


% plotting

t_1_2 = vertcat(t1,t2);
y_1_2 = vertcat(y1,y2);

subplot(2,3,1)
plot(t_1_2, y_1_2(:,2)); 
xlabel('t')
ylabel('g/L')
legend('S')

subplot(2,3,2)
plot(t_1_2, 1000*y_1_2(:,[1,3])); 
xlabel('t')
ylabel('mg/L')
legend('X','P')

subplot(2,3,3)
plot(t_1_2, y_1_2(:,4:5)); 
xlabel('t')
legend('Cx','Nx')

subplot(2,3,4)
plot(t_1_2, y_1_2(:,6:10)); 
xlabel('t')
ylabel('mmol/L')
legend('NT','QT','RT','PT','ST')

subplot(2,3,5)
plot(t_1_2, y_1_2(:,11:15)); 
xlabel('t')
ylabel('mmol/L')
legend('K+','Mg2+','Ca2+','M2+','Cl-')

subplot(2,3,6)
plot(t_1_2, y_1_2(:,16)); 
xlabel('t')
ylabel('mmol/L')
legend('MAP')












