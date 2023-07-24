% medium components

% pH and precipitation


function F = pH_ver_1(x)
% x(1) = OH-            % x(8) = H2SO4          
% x(2) = PO4(3-)        % x(9) = gamma_-3       
% x(3) = HPO4(2-)       % x(10) = gamma_-2      
% x(4) = H2PO4-         % x(11) = gamma_-1    
% x(5) = H3PO4          % x(12) = gamma_0
% x(6) = SO4(2-)        % x(13) = gamma_1
% x(7) = HSO4-          % x(14) = I

H = 10^(-6.5);
A = 0.51;
y = zeros(16,1);
y(6) = 1;
y(12) = 1;
y(11) = 1;
y(9) = 1;
y(10) = 1;

Ka3_1 = 1;
Ka2_1 = 2;
Ka1_1 = 3;
Ka2_2 = 1;
Ka1_2 = 2.1;

% charge balance
F(1) = H-x(1)-3*x(2)-2*x(3)-x(4)-2*x(6)-x(7)+y(6)+2*y(12)+y(11);

% acid-base equilibrium
F(2) = -Ka3_1 + (H*x(2)/x(3))*(x(13)*x(9)/x(10));
F(3) = -Ka2_1 + (H*x(3)/x(4))*(x(13)*x(10)/x(11));
F(4) = -Ka1_1 + (H*x(4)/x(5))*(x(13)*x(11)/x(12));
F(5) = -Ka2_2 + (H*x(2)/x(3))*(x(13)*x(10)/x(11));
F(6) = -Ka1_2 + (H*x(2)/x(3))*(x(13)*x(11)/x(12));

% davis equation
F(7) = -x(9) + 10^(-A*9*(sqrt(x(14))/(1+sqrt(x(14))) - 0.3*x(14)));
F(8) = -x(10) + 10^(-A*4*(sqrt(x(14))/(1+sqrt(x(14))) - 0.3*x(14)));
F(9) = -x(11) + 10^(-A*1*(sqrt(x(14))/(1+sqrt(x(14))) - 0.3*x(14)));
F(10) = -x(12) + 1;
F(11) = -x(13) + 10^(-A*1*(sqrt(x(14))/(1+sqrt(x(14))) - 0.3*x(14)));
F(12) = -x(14) + 0.5*(9*x(2)+4*x(3)+x(4)+4*x(6)+x(7)+y(6)+4*y(12)+y(11));

% total phosphate, sulfate species equation
F(13) = y(9) + x(2)+x(3)+x(4)+x(5);
F(14) = y(10) + x(6)+x(7)+x(8);