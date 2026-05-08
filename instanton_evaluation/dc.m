% O(h^4) center-difference derivative
% f' = (-f(t+2h) + 8f(t+h) - 8f(t-h) + f(t-2h))/(12h)

function grad = dc(x,h)

grad = (-x(5:end)+8*x(4:end-1)-8*x(2:end-3)+x(1:end-4))/(12*h);
grad1 = (-3*x(1)+4*x(2)-x(3))/(2*h);
grad2 = (-3*x(2)+4*x(3)-x(4))/(2*h);
gradend_1 = (3*x(end-1)-4*x(end-2)+x(end-3))/(2*h);
gradend =  (3*x(end)-4*x(end-1)+x(end-2))/(2*h);
% dc1= interp1([1:5],[grad2 grad(1:4)],0,'spline');
% dcend= interp1(1:5,[grad(end-3:end) gradT_1],6,'spline');

grad=[grad1, grad2, grad, gradend_1, gradend];

end

