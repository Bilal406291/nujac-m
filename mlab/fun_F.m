
function F = fun_F(p,N,y_old,dt,ep,tau)

 F = (p-y_old) - dt*f(p,N,ep,tau);
 
end