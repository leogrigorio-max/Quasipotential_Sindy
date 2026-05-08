function psymp = int_symplectic_p(C,E0,x,p,dt)

dt = -dt;  %backwards integration
N = length(x);
for n=(N-1):-1:1
    
%     v = INR_p_v1([x(:,n);p(:,n)],[x(:,n+1);p(:,n+1)],dt);
    v = INR_p(C,E0,x(1,n),x(2,n),p(1,n),p(2,n),x(1,n+1),x(2,n+1),p(1,n+1),p(2,n+1),dt);
    p(1,n) = v(1); p(2,n) = v(2);
    
end
    
psymp = p;    
    
end
