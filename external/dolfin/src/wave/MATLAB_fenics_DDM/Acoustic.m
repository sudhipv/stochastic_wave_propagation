function integral = Acoustic(u, v, w, c, du, dv, dw, dx, ds, x, d, t, eq)

% Variational formulation of Poisson's equation.
% Copyright (C) 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.

if eq == 1
  integral = c*c*du'*dv*dx + g(x,d,t)*u*v*ds;
  
elseif eq == 2
  integral =  u*v*dx + g(x,d,t)*u*v*ds;
  
else       
  integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;
end 

%--- Conductivity (penalty factor) ---
function y = g(x, d, t)
%y = 1000;
y = 1e7;

%--- Dirichlet boundary condition ----
function y = gd(x, d, t)
y = 0;

%--- Neumann boundary condition ---
function y = gn(x, d, t)
y = 0;

%--- Right-hand side, source term ---
function y = f(x, d, t)
   y = 0;
%  y = 100*exp(-0.01*t*((x(1)-0.5).^2+(x(2)-0.5).^2));
% y = 40*exp(-100*(x(1)-0.5).^2+(x(2)-0.5).^2);
% y = 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
% y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));
