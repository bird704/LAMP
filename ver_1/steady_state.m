% MODEL DEVELOPMENT

% Use 4th order Runge-Kutta method
% to solve the systems of 1st order - IV ODE problem


function c_vector = steady_state(h, L, D, U, k, c_in)
n = (L-0)/h;   % number of intervals
e = ones(n+1,1); % number of points = n+1

alpha = (U*h)/(2*D);
beta = (k*h^2)/(2*D);
a = e.*alpha; b = e.*beta;

% build the tridiagonal matrix
A = spdiags([-e-a, 2.*(e+b), -e+a],[-1 0 1],n+1,n+1);
A(1,1) = 4*alpha*(1+alpha)+2*(1+beta); A(1,2) = -2;
A(n+1,n) = -2; A(n+1,n+1) = 2*(1+beta);

% build the f vector
f=zeros(n+1,1); f(1,1) = 4*alpha*(1+alpha)*c_in;

% solve the linear system
c_vector = A\f;

