% ==  Stepwise Sparse Regressor (SSR) Sindy algorithm  applied to instanton data  ==
% F => data matrix we want to find functional form: F = gradient or F = transverse drift component
% F = Theta*Xi, Theta contains the dictionary matrix

% sindy method is then applied to infer the underlying structure of drift
% b(x) = -grad V(x) + l(x); p = 2 grad V(x)
%  From sindy we obtain how grad V(x) depends on x
% also we obtain the transverse component l(x) = \dot \phi(t) - p(t)/2



% == reading data (instanton) ==
file = 'trajectories/limit_cycle_lambda_-15.000.mat';

load(file)
 Z = [x(1,:); x(2,:); p(1,:); p(2,:); dotX(1,:); dotX(2,:)];

toler = 1e-2;

interv = 1:length(p(1,:));

K = 6; %number of partition sets - one for training, the rest for control
siz_int = floor(interv(end)/K);


% ==== chosing F=-p/2 to learn the gradient or F = dotX' - p'/2; to learn the orthogonal component
F = -p'/2;  %learning gradient:  recall that  p = 2 grad V
% F = dotX' - p'/2; % learning transverse component field


% == building library ==
% for paper example 4.1 we recommend restricting the library to monomials up to 3rd order to get sparse solution
library = @(x,y,p,q,xdot,ydot) [ones(length(x),1), x, y, x.^2, y.^2, x.*y, x.^3, y.^3, x.^2.*y, x.*y.^2, x.^4, y.^4, x.^3.*y, x.^2.*y.^2, x.*y.^3];
lib = "1, x, y, x.^2, y.^2, x.*y, x.^3, y.^3, x.^2.*y, x.*y.^2, x.^4, y.^4, x.^3.*y, x.^2.*y.^2, x.*y.^3";
lib = split(lib,", ");

% == score: quadratic error between training and validation
score = zeros(K,2*length(lib)); %2=columns of F

% == K-fold cross validation
for n = 1:K

interv_c = interv(interv<=(n*siz_int) & (n-1)*siz_int<interv); %indices for control
interv_t = interv(not(interv<=(n*siz_int) & (n-1)*siz_int<interv)); %indices for training

Ft = F(interv_t,:); %training
Fc = F(interv_c,:); %control

[m,c] = size(Ft); [mc,cc] = size(Fc);

Theta = library(Z(1,interv_t)',Z(2,interv_t)',Z(3,interv_t)',Z(4,interv_t)',Z(5,interv_t)',Z(6,interv_t)'); %training library
Theta_c = library(Z(1,interv_c)',Z(2,interv_c)',Z(3,interv_c)',Z(4,interv_c)',Z(5,interv_t)',Z(6,interv_t)'); %control library

[mTh,cTh] = size(Theta);

% === Stepwise Sparse Regressor (SSR) ===

Xi_1st = Theta\Ft; % initial guess: least-squares
Xi = Xi_1st;
Xir = reshape(Xi,1,c*cTh);
[fv, id] = min(abs(Xir));

Xis = zeros(cTh,c,c*cTh); %set of Xi's with different sparsity levels
for q=1:c*cTh  % q: level of sparsity - number of zero coeffs
  % find minimum coefficient and sets it to zero

  if (mod(id,cTh)==0)
    id_i=cTh;
  else
    id_i=mod(id,cTh);
  end
  id_j = floor((id-1)/cTh)+1; %argmin matrix indices
  Xi(id_i,id_j) = 0;
  zeroinds = (abs(Xi)<toler);
  for col = 1:c
    non_zinds = ~zeroinds(:,col);
    Xi(non_zinds,col) = Theta(:,non_zinds)\Ft(:,col);
  end
  Xir = reshape(Xi,1,c*cTh);
  Xi_nozero = Xir(~Xir==0);
  if ~isempty(Xi_nozero)
    fv = min(abs(Xi_nozero));
    id = find(abs(Xir)==fv);
  end
  % id
  % Xi
  Xis(:,:,q) = Xi;
  score(n,q) = trace((Fc-Theta_c*Xi)'*(Fc-Theta_c*Xi)); %testing against control data

end

end


final_score = mean(score,1);
%  == (balance between sparsity and acuracy || under/over-fitting)
[sc, j] = min(final_score);
best_Xi = Xis(:,:,j);  %if best_Xi is not sparse, look for Pareto front - see our appendix to get the sparsest possible model

hold off
semilogy(1:q,final_score,'.','markersize',20)
xlabel('sparsity $q$','interpreter','latex')
title( 'Cross validation score','interpreter','latex')
set(gca,'GridLineStyle','--')
grid on

% === printing results: active terms obtained by sparse learning ===
fprintf('**** sparse learning results ***** \n')
fprintf('\n  library   ||   p1 coeffs   ||   p2 coeffs \n')
disp([lib,best_Xi])
