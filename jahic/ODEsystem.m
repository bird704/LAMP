% model development

% important variables : 
% X = concentration of biomass
% S = concentration of substrate
% P = concentration of protein
% Cx = concentration of intracellular(X) carbon
% Nx = concentration of intracellular(X) nitrogen

function f = ODEsystem(t, y, p, M)

%% celluar metabolism
% set parameters for cellular metabolism
F_in = 0.5;                 % p(1) = qm
V_L = 1;                    % p(2) = Yem
y1_in = 0;                  % p(3) = qs_max
y2_in = 1;                  % p(4) = Ks
Cp = 1;                     % p(5) = KaE
Np = 1;                     % p(6) = H       
                            % p(7) = a_P
                            % p(8) = a_C
                            % p(9) = a_N

% explicit equations                   
qs = p(3)*( y(2) / ( p(4) + y(2) ) );
mu = ( qs - p(1) )*p(2);
qp = ( p(7) * mu ) / ( 10^( -log10(p(5)) + log10(p(6)) ) + 1 );
qcx = p(8)*mu;
qnx = p(9)*mu;

% differential equations for cellular metabolism
f = zeros(5,1);
f(1) = (F_in/V_L)*( y1_in - y(1) ) + mu*y(1);   % y(1) = X
f(2) = (F_in/V_L)*( y2_in - y(2) ) - qs*y(1);   % y(2) = S
f(3) = -(F_in/V_L)*y(3) + qp*y(1);              % y(3) = P
f(4) = qcx - mu*y(4) - qp*Cp;                   % y(4) = Cx
f(5) = qnx - mu*y(5) - qp*Np;                   % y(5) = Nx

%% medium components
% set parameters for medium components
A = F_in/V_L;
M_in = [1,1,1,1,1,1,1,1,1,1];
qM = [1,1,1,1,1,1,1,1,1,1];
rpM = [1,1,1,1,1,1,1,1,1,1];
Mx = [1,1,1,1,1,1,1];
nM = 5;

% explicit equations
for i = 1:3
    qM(i) = (nM*y(i+5)/(y(6)+2*y(7)+4*y(8)))*qnx;
end
for i = 4:10
    qM(i) = Mx(i-3)*mu;
end

% differential equations for medium components
f(6) = A*(M_in(1)-y(6)) - qM(1)*y(1) - rpM(1);      % y(6)= ammonium(NT)
f(7) = A*(M_in(2)-y(7)) - qM(2)*y(1) - rpM(2);      % y(7)= glutamine(QT)
f(8) = A*(M_in(3)-y(8)) - qM(3)*y(1) - rpM(3);      % y(8)= arginine(RT)
f(9) = A*(M_in(4)-y(9)) - qM(4)*y(1) - rpM(4);      % y(9)= phosphate(PT)
f(10) = A*(M_in(5)-y(10)) - qM(5)*y(1) - rpM(5);    % y(10)= sulfate(ST)
f(11) = A*(M_in(6)-y(11)) - qM(6)*y(1) - rpM(6);    % y(11)= potassium(+)
f(12) = A*(M_in(7)-y(12)) - qM(7)*y(1) - rpM(7);    % y(12)= magnesium(2+)
f(13) = A*(M_in(8)-y(13)) - qM(8)*y(1) - rpM(8);    % y(13)= calcium(2+)
f(14) = A*(M_in(9)-y(14)) - qM(9)*y(1) - rpM(9);    % y(14)= trans M(2+)
f(15) = A*(M_in(10)-y(15)) - qM(10)*y(1) - rpM(10); % y(15)= chloride(-)


end