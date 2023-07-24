% medium components

% pH and precipitation


function F = pH_ver_1(x)
% x(1) = OH-            % x(8) = H2SO4          
% x(2) = PO4(3-)              
% x(3) = HPO4(2-)        
% x(4) = H2PO4-          
% x(5) = H3PO4          
% x(6) = SO4(2-)        
% x(7) = HSO4-          

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
F(2) = -Ka3_1 + (H*x(2)/x(3));
F(3) = -Ka2_1 + (H*x(3)/x(4));
F(4) = -Ka1_1 + (H*x(4)/x(5));
F(5) = -Ka2_2 + (H*x(2)/x(3));
F(6) = -Ka1_2 + (H*x(2)/x(3));

% total phosphate, sulfate species equation
F(13) = y(9) + x(2)+x(3)+x(4)+x(5);
F(14) = y(10) + x(6)+x(7)+x(8);