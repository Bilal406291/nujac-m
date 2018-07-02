function J=jacNum(g,x)

% computes the Jacobian of a function

n=length(x);

J=zeros(n,n);

%fx=feval(func,x);

fx=g(x);

ep=1.e-8;  % could be made better

xx=x;

for i=1:n
    
xx(i)=xx(i)+ep;

J(:,i)=(g(xx)-fx)/ep;

%J(:,i)=(feval(func,xx)-fx)/eps;

xx(i)=x(i);

end