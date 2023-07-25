% model development

% important variables : 
% X = concentration of biomass
% S = concentration of substrate
% OCR = oxygen consumption rate

function f = ODEsystem_jahic(t, y, p, split, muu)
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
         
% Feed_rate (mL/h)
F_in = feed_rate(t,split, muu);

% Volume (L)
V = 3;

% explicit equations                   
qs = p(5) * ( y(2) / ( p(3) + y(2) ) );
qs_an = ( qs - p(4) ) * p(7) * ( p(2) / p(1) );
qs_en = ( qs - qs_an );
qo = qs_an * p(8) + qs_en * p(9);
mu = ( qs - p(4) ) * p(7);

% DAEs for cellular metabolism
f = zeros(3,1);
f(1) = ( -F_in / V + mu ) * y(1);                     % y(1) = X
f(2) = ( F_in / V ) * ( p(6) - y(2) ) - qs * y(1);    % y(2) = S
f(3) = - y(3) + qo * y(1) * V;                         % y(3) = OCR
end