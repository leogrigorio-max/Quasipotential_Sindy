%================= Instanton Newton-Raphson  ===================
% Newton-Raphson method to solve implicit ODE algebraic equation
% Both function F and Jacobian should be specified in each case
% F is the function to be minimized (zero) whose roots we want
% see symplectic.m file for jacobian
%========================  ===================

function NR = INR_x(C,E0,x,y,px,py,x_,y_,px_,py_,h)
    u = [x;y];
    NR = [0;0];
    tol = 1.e-9; %tolerance (for convergence)


    for j = 1:800

     %    Jacobian
     J = jacobian_x_limit_cycle(C,E0,h,px,px_,py,py_,x,x_,y,y_);
     f = hamilton_x_limit_cycle(C,E0,h,px,px_,py,py_,x,x_,y,y_);
    
     u = u - J\f;
     x = u(1); y = u(2);

     delta = sqrt(f'*f);
%        disp(delta)
        if (delta<tol)
            NR(1) = u(1); NR(2) = u(2);
%            disp(f)
             % disp(j)
            break
        end
    end
if (j==800)
    sprintf('INR x did not converge')
end

end
