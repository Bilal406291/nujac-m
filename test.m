

N=19;

h = 1/(N+1);

xint = h*(1:N)';

uint = sin(2*pi*xint) + (1/2)*sin(pi*xint);

u0=[uint; xint];

%uv=Rho_test(u0,N);

M=mass(u0);

rhs=f(u0);

%myf = @(x) f(x);

jac1=jacNum(@f,u0);

%jac2=NumJacob(@f,u0);

%jac3=jacobFD(@f,u0);

q=M\rhs;




