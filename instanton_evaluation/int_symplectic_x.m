function xsymp = int_symplectic_x(C,E0,x,p,dt)

N = length(x);
for n=1:(N-1)
    
%     v = INR_x_v1([x(:,n+1); p(:,n+1)],[x(:,n);p(:,n)],dt);
    v = INR_x(C,E0,x(1,n+1),x(2,n+1),p(1,n+1),p(2,n+1),x(1,n),x(2,n),p(1,n),p(2,n),dt);
    x(1,n+1) = v(1); x(2,n+1) = v(2);
    
end
    
xsymp = x;    
    
end
