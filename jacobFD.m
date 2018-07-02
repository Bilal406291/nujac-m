 function J = jacobFD(g,x) 
% Calculates the Jacobian of the
% system of non-linear equations:
% f(x) = 0, through finite differences.
% The Jacobian is built by columns

delx=1e-8;

m=length(x);

J=zeros(m,m);

for j = 1:m 
   xx = x; 
   xx(j) = x(j) + delx; 
   J(:,j) = (g(xx)-g(x))/delx; 
end
