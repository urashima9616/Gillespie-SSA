function [A,t] = gillespiessa_anim(Nanim, knu, c,p, A0,tfin)
% plots the results after each Nanim steps
% Gillespe Stochastic Simulation Algorithm for a system of q chemical reactions
% Nchs ... number of chemical species
% q ... number of chemical reactions
%
% i-th reaction with the rate constant knu(i), i=1,2,...,q:
%   c(1,i) A(1) + c(2,i) A(2) + ... + c(Nchs,q) A(Nchs)
%    --->
%   p(1,i) A(1) + p(2,i) A(2) + ... + p(Nchs,q) A(Nchs)
%
%  INPUT:
%   Nanim ... plots the results after each Nanim steps
%    knu ... 1-by-q row vector; rate constants of the chemical reactions;
%            knu = k/nu^(order-1), where nu=volume and order=order of the chemical reaction
%    c ... Nchs-by-q; c(j,i) = coeficient at the j-th chemical species at the reactant side of the i-th reaction
%    p ... Nchs-by-q; p(j,i) = coeficient at the j-th chemical species at the product side of the i-th reaction
%    A0 ... Nchs-by-1 column vector;
%           A0(j) = number of molecules of j-th chemical species at the initial time t=0
%    tfin ... final time; if the simulation raches tfin it stops
%
%  OUTPUT:
%    A ... Nchs-by-N; A(j,m) = number of j-molecules in time interval [t(m), t(m+1)]
%          A(:,1) == A0(:);
%    t ... 1-by-N;


% Gillespie SSA
% =============
% (a) generate two randon numbers r1 and r2 uniformly distributed in (0,1)
% (b) Compute alpha0 = \sum_{i=1}^q alpha_i(t),
%     where alpha_i(t) = propensity function of i-th reaction
% (c) The next reaction take place at time t+tau, where
%     tau = 1/alpha0 * log(1/r1)   ... natural logrithm
% (d) Determine which reaction occurs. Find j such that
%     r2 >= 1/alpha0 \sum_{i=1}^{j-1} alpha_i(t)
%     r2 <  1/alpha0 \sum_{i=1}^j     alpha_i(t)
%     Update the number of reactants and products of the j-th reaction
% GO TO (a) with t = t + tau

% Initial checks
if (size(A0,2) > 1)
  error('The initial numbers of molecules A0 must form a COLUMN vector');
end
Nchs = size(A0,1);

if (size(knu,1) > 1)
  error('The rate constants knu must form a ROW vector');
end
q = size(knu,2);

if (size(c) ~= [Nchs,q])
  error('The reactant coeficients must be Nchs-by-q but Nchs=%d, q=%d and size(c)=(%d,%d)', Nchs,q,size(c,1),size(c,2));
end

if (size(p) ~= [Nchs,q])
  error('The product coeficients must be Nchs-by-q but Nchs=%d, q=%d and size(p)=(%d,%d)', Nchs,q,size(p,1),size(p,2));
end

% assume maximum of Nalloc steps of the algorithm
% allocate memory for Nalloc steps
Nalloc = 10000;

A = zeros(Nchs,Nalloc);
t = zeros(1,Nalloc);

% set the initial values
m = 1;
A(:,1) = A0;
t(1) = 0;

% zero-order chemical reaction:
%    0 -> A; alpha = k\nu;  A(t) = A(t) + 1
% first-order chemical reaction:
%    A -> 0; alpha = A(t)*k;  A(t) = A(t) - 1

% Propensity functions of chemical reactions:

%  order = zeros(1,q);
%  for j = 1 : q
%    order(j) = sum(c(:,j));
%    if (order(j) == 0)
%    end % if
%  end % for j

% initialize the random number generator
rand('twister',sum(1000000*clock));

% Gillespie SSA
while (t(m) <= tfin)
  % (a)
  r = rand(1,2); % 1-by-2 vector of two random numbers
  % (b)
  alpha = zeros(1,q);
  for i=1:q % go thru all reactions and compute the propensity function alpha(i)
    alpha(i) = knu(i);
    for j=1:Nchs % go thru all chemical species
      for s=1:c(j,i)  % ... alphai=knu*A*(A-1)*(A-2)*...*(A-c(j,i)+1)
        alpha(i) = alpha(i) * (A(j,m)-s+1);
      end % for s
    end % for j
    %alpha0 = alpha0 + alpha(i);
  end % for i
  alpha0 = sum(alpha);

% -------- CHECK for test2 --------
%  TOL=1e-12;
%  alpha0_test = A(1,m)*knu(1) + knu(2);
%  if (abs(alpha0-alpha0_test)>TOL) error('test1 failed a0=%g  a0_test=%g',alpha0,alpha0_test); end
% -------- END CHECK --------------

% -------- CHECK for test3 --------
%  TOL=1e-12;
%  alpha1_test = A(1,m)*(A(1,m)-1)*(A(1,m)-2)*knu(1);
%  alpha2_test = A(1,m)*(A(1,m)-1)*knu(2);
%  alpha3_test = A(1,m)*knu(3);
%  alpha4_test = knu(4);
%  alpha0_test = alpha1_test + alpha2_test + alpha3_test + alpha4_test;
%  if (abs(alpha(1)-alpha1_test)>TOL*alpha1_test)
%    error('test1 failed a1=%g  a1_test=%g  abserr=%g',alpha(1),alpha1_test, abs(alpha(1)-alpha1_test));
%  end
%  if (abs(alpha(2)-alpha2_test)>TOL*alpha2_test)
%    error('test1 failed a2=%g  a2_test=%g  abserr=%g',alpha(2),alpha2_test, abs(alpha(2)-alpha2_test));
%  end
%  if (abs(alpha(3)-alpha3_test)>TOL*alpha3_test)
%    error('test1 failed a3=%g  a3_test=%g  abserr=%g',alpha(3),alpha3_test, abs(alpha(3)-alpha3_test));
%  end
%  if (abs(alpha(4)-alpha4_test)>TOL*alpha4_test)
%    error('test1 failed a4=%g  a4_test=%g  abserr=%g',alpha(4),alpha4_test, abs(alpha(4)-alpha4_test));
%  end
%  if (abs(alpha0-alpha0_test)>TOL*alpha0)
%    error('test1 failed a0=%g  a0_test=%g  abserr=%g', alpha0, alpha0_test, abs(alpha0-alpha0_test));
%  end
% -------- END CHECK --------------

  alpha0_1 = 1/alpha0;
  % (c)
  tau = alpha0_1 * log(1/r(1));
  %(d)
  i = 1;
  sumalpha = alpha0_1*alpha(1);  % == 1/alpha0 \sum_{i=1}^{j} alpha_i(t)
  while (r(2) >= sumalpha)
    i = i + 1;
    sumalpha = sumalpha + alpha0_1*alpha(i);
  end % while
  % i-th reaction occurs
  %fprintf('m=%d A(m)=%d t(m)=%f i=%d',m,A(1,m),t(m),i);
  A(:,m+1) = A(:,m) - c(:,i) + p(:,i);
  t(m+1) = t(m) + tau;
  %fprintf(' t(m+1)=%f A(m+1)=%d  alpha=[%g,%g,%g,%g]\n',t(m+1),A(1,m+1),alpha(1),alpha(2),alpha(3),alpha(4));

% -------- CHECK for test2 --------
%  if (r(2) < A(1,m)*knu(1)/alpha0)
%    Atest = A(1,m) - 1;
%  else
%    Atest = A(1,m) + 1;
%  end
%  if (abs(A(1,m+1)-Atest)>TOL)
%    error('test2 failed A=%d Atest=%d',A(1,m+1),Atest);
%  end
% -------- END CHECK --------------

%  % -------- CHECK for test3 --------
%  if (r(2) < alpha1_test/alpha0)
%    Atest = A(1,m) - 1;  i_test=1;
%  elseif (r(2) < (alpha1_test+alpha2_test)/alpha0)
%    Atest = A(1,m) + 1; i_test=2;
%  elseif (r(2) < (alpha1_test+alpha2_test+alpha3_test)/alpha0)
%    Atest = A(1,m) - 1; i_test=3;
%  else
%    Atest = A(1,m) + 1; i_test=4;
%  end
%  if (i ~= i_test)
%    error('test2 failed i=%d i_test=%d',i,i_test);
%  end
%  if (abs(A(1,m+1)-Atest)>TOL)
%    error('test2 failed A=%d Atest=%d',A(1,m+1),Atest);
%  end
%  % -------- END CHECK --------------


  % Animation
  if (mod(m,Nanim) == 0)
    %figure(1);
    subplot(1,2,1);
    plot(t(1:m), A(1,1:m), t(m),A(1,m),'.r');
    xlabel('time');
    ylabel('number of X molecules');
    axis( [0 tfin 0 700]);

    %figure(2);
    subplot(1,2,2);
    plot(A(1,1:m),A(2,1:m), A(1,m),A(2,m),'.r');
    xlabel('number of X molecules');
    ylabel('number of Y molecules');
    axis( [0 700 200 1400]);
    drawnow;
  end

  m = m + 1;

  % memory management
  if (m >= Nalloc)
    fprintf('Doubling memory for data. (Reallocation.)  t(%d)=%f', m, t(m));
    tic
    Aaux = zeros(Nchs,2*Nalloc);
    taux = zeros(1,2*Nalloc);
    Aaux(:,1:Nalloc) = A;
    taux(1:Nalloc) = t;
    A = Aaux;
    t = taux;
    clear Aaux taux;
    Nalloc = 2*Nalloc;
    fprintf('  done. [%fsec]\n', toc);
  end
end % while


% cutting the zeros at the end of arrays
A = A(:,1:m);
t = t(1:m);

% postprocessing
if (t(m) > tfin)
  t(m) = tfin;
end
for j=1:Nchs
  if (A(j,m) < 0)
    A(j,m) = 0;
  end
end

