function H = hamiltonian(x,p,N)

H = zeros(1,N);

% B = @(x,y) [-x.^3+y.^3;-y.^3-x.^3];
% B = @(x,y) [-x.*(x.^2+y.^2-1) + y; -y.*(x.^2+y.^2-1) - x];
B = @(x,y) [-x+y;-y-x];
a = @(x,y) eye(2)./(1+x^2);
for n=1:N

    H(n) = 0.5*p(:,n)'*a(x(1,n),x(2,n))*p(:,n) + p(:,n)'*B(x(1,n),x(2,n));

end
% H = (1/2)*<p,Dp> + <b,p>; rescaled, D = eps Id
% eps H =  \bar H = 1/2 <p,p> + <b,p>
end
