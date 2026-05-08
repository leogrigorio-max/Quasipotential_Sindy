% % =====================================================
% Instanton calculation
% % =====================================================

% tic
itmax = 400; % maximum number of iterations
alpha0 = .1;  %sometimes should be small as 0.06
tol = 1.e-9;   %tolerance
E0 = 1.e-7;

% == instanton initial condition - starting at attractor (limit cycle - radius = 1)
x0 = [-1.0; 0.0];
x0 = x0/norm(x0);
x(:,1) = x0;
x_(:,1)= x0;

%======== Starting Iterations ======

lambda = -12;   %final condition for p; should be increased incremently 

for it=1:itmax

  z = x(1,end);

  p(1,end) = lambda;
% if (it<20)
%    p(1,end) = lambda;
% else
% p(1,end) = lambda/(z);  % F = +-log|z|
  % p(1,end) = lambda/(z*log(z));  % F = +-loglog|z|
% end

    p_ = p;
    x_ = x;
    if (it==itmax)
		disp('Number max of iteractions reached!')
		break
    end

	xold = x(1,end);

   % =====  Integrating Eqns of motion ====

   p = int_symplectic_p(C,E0,x,p,ds);
   x = int_symplectic_x(C,E0,x,p,ds);

   x = alpha0*x + (1.0 - alpha0)*x_;
   p = alpha0*p + (1.0 - alpha0)*p_;


   % ==== Computing arclength =====
   dx1 = dc(x(1,:),ds); dx2 = dc(x(2,:),ds);
   mod_dx = sqrt(dx1.*dx1 + dx2.*dx2);
   C = trapz(s,mod_dx)/s(end);
   % =============================

  E = hamiltonian(x,p,N);

  if (((isnan(x(1,end)))||abs(xold)<tol)&&(it>8))
	   disp('Nan or vanishing solution')
	   break
	end
	delta = abs(x(1,end)-xold)/abs(xold);
  if((delta < tol) && (it>2))
      break
  end

% ============= Monitoring error ============
monitor(it,x(1,end),E,mod_dx,delta,s);

end  % loop it


% elapse = toc
% ======= Plotting instanton =====
figure(3)
set(gcf,'position',[650 120 550 500])
plot(x(1,:),x(2,:),'-.')

%  === saving data ==
instanton_file = sprintf('trajectories/limit_cycle_lambda_%.3f.mat',lambda)
% save(instanton_file)
