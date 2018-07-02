

clear; clc;

format long

ti=0.0;

tf=1.0;

nt=100;

dt=(tf-ti)/nt;

ep=1e-2;

tau = 1e-2;

N=19;

h = 1/(N+1);

xint = h*(1:N)';

uint = sin(2*pi*xint) + (1/2)*sin(pi*xint);

u0=[uint; xint];

Y=zeros(2*N,nt);

Y(:,1)=u0;

max_iter=100;

tol=1e-8;

 v=Y(:,1);
 
 k=zeros(nt,1);
 
for m=2:nt
    
    it=0; % counter
    
    t=ti+m*dt;
    
    % newton loop 
  
    for p=1:max_iter 
        
        Fn=fun_F(v,N,Y(:,m-1),dt,ep,tau);
        
        Jn=jacobFD(@fun_F,v,N,Y(:,m-1),dt,ep,tau);
        
        ls=Jn\Fn;
        
        z=v-ls;
        
        it=it+1;
       
        if norm(abs(z-v))<tol 
             
             disp('converged')
             
             break
             
         end
         
        v=z;
        
    end
    
    k(m)=it;
    
    Y(:,m) = z;
end  

Y(:,end)   






