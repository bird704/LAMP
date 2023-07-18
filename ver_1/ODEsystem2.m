% model development

% important variables : 
% X = concentration of biomass
% S = concentration of substrate
% P = concentration of protein
% Cx = concentration of intracellular(X) carbon
% Nx = concentration of intracellular(X) nitrogen

function f = ODEsystem2(t, y, p, k, M)
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
         
% F_in
if t < 20
    F_in = 0.0385;
else
    F_in = 0.0105;
end


% explicit equations                   
qs = p(3)*( y(2) / ( p(4) + y(2) ) );
mu = ( qs - p(1) )*p(2);
qp = ( p(7) * mu ) / ( 10^( -log10(p(5)) + log10(p(6)) ) + 1 );
qcx = p(8)*mu;
qnx = p(9)*mu;

test = -(F_in/k(1))+mu

% differential equations for cellular metabolism
f = zeros(5,1);
f(1) = (F_in/k(1))*( k(2) - y(1) ) + mu*y(1);   % y(1) = X
f(2) = (F_in/k(1))*( k(3) - y(2) ) - qs*y(1);   % y(2) = S
f(3) = -(F_in/k(1))*y(3) + qp*y(1);              % y(3) = P
f(4) = qcx - mu*y(4) - qp*k(4);                   % y(4) = Cx
f(5) = qnx - mu*y(5) - qp*k(5);                   % y(5) = Nx



%% medium components
% set parameters for medium components
A = F_in/k(1);
M_in = [1,1,1,1,1,1,1,1,1,1,0];
qM = [1,1,1,1,1,1,1,1,1,1];
rpM = zeros(11,1);
Mx = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
nM = [1,2,4];

% explicit equations
for i = 1:3
    qM(i) = (nM(i)*y(i+5)/(y(6)+2*y(7)+4*y(8)))*qnx;
end
for i = 4:10
    qM(i) = Mx(i-3)*mu;
end

k1 = 1.3*10^(11);
Ksp1 = 10^(-12.24);
Qsp1 = y(12)*y(6)*y(9);

if Qsp1 < Ksp1
    rpM(1,4,7,11) = 0;
else
    rpM(1,4,7,11) = k1*(Qsp1 - Ksp1);
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
f(16) = A*(M_in(11)-y(16)) + rpM(11);               % y(16)= MAP

%% solving pH model
fun = @pH_ver_1;
x0 = [1,2,11,13,11,2,1,4,1,15,1,1,16,1];
%x = fsolve(fun,x0)

end