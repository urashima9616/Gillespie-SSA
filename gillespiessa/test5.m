% Chemical system with SNIPER bifurcation:
%       0 --> Y      with rate constant k1d
%       Y --> X      with rate constant k2d
%       X --> 0      with rate constant k3d
%      2X --> 3X     with rate constant k4d
%      3X --> 2X     with rate constant k5d
%   X + Y --> X + 2Y with rate constant k6d
%  2X + Y --> 2X     with rate constant k7d

% Nchs ... number of chemical species
% q ... number of chemical reactions
%    knu ... 1-by-q row vector; rate constants of the chemical reactions;
%            knu = k/nu^(order-1), where nu=volume and order=order of the chemical reaction
%    c ... Nchs-by-q; c(j,i) = coeficient at the j-th chemical species at the reactant side of the i-th reaction
%    p ... Nchs-by-q; p(j,i) = coeficient at the j-th chemical species at the product side of the i-th reaction
%    A0 ... Nchs-by-1 column vector;
%           A0(j) = number of molecules of j-th chemical species at the initial time t=0
%    tfin ... final time; if the simulation raches tfin it stops

%clear all;

k1d = 12;   % sec^(-1) mm^(-3)
k2d = 1;    % sec^(-1)
k3d = 33;   % sec^(-1)
k4d = 11;   % sec^(-1) mm^3
k5d = 1;    % sec^(-1) mm^6
k6d = 0.6;  % sec^(-1) mm^3
k7d = 0.13; % sec^(-1) mm^6
vol = 40; % mm^3
k1nu = k1d*vol;   % sec^(-1)
k2nu = k2d;       % sec^(-1)
k3nu = k3d;       % sec^(-1)
k4nu = k4d/vol;   % sec^(-1)
k5nu = k5d/vol^2; % sec^(-1)
k6nu = k6d/vol;   % sec^(-1)
k7nu = k7d/vol^2; % sec^(-1)

q = 7;
Nchs = 2;

knu = [ k1nu k2nu k3nu k4nu k5nu k6nu k7nu ];
c = [0, 0, 1, 2, 3, 1, 2;
     0, 1, 0, 0, 0, 1, 1];
p = [0, 1, 0, 3, 2, 1, 2;
     1, 0, 0, 0, 0, 2, 0];
A0 = [70;
      500];
tfin = 100;

%ode45

% plot Nreal realizations
Nreal = 1;
figure(1)
clf;
%hold on;
figure(2)
clf;
%hold on;
for n = 1 : Nreal
  [A,t] = gillespiessa(knu, c,p, A0,tfin);

  %save(['test5_A_r' num2str(Nsave+n) '.dat'], 'A', '-ascii');
  %save(['test5_t_r' num2str(Nsave+n) '.dat'], 't', '-ascii');
  save(['test5_A.dat'], 'A', '-ascii');
  save(['test5_t.dat'], 't', '-ascii');

  fprintf('Plotting...');
  [Apl,tpl] = mkplotdata(A,t);
  Ac{n} = Apl;
  tc{n} = tpl;
  co = get(gca,'ColorOrder');
  figure(1)
  plot( tpl,Apl(1,:), 'Color', co(mod(n,7)+1,:));
  figure(2)
  plot( A(1,:),A(2,:), 'Color', co(mod(n,7)+1,:));

  fprintf('done [%fsec].\n', toc);
end % for n

%  tm = [0:0.1:tfin];
%  Mm = k2nu/k1*(1-exp(-k1.*tm));
%  plot( tm,Mm, 'k--', 'LineWidth',3);

figure(1)
%hold off;
axis( [0 tfin 0 400]);
xlabel('time');
ylabel('number of X molecules');
title('Chemical system with SNIPER bifurcation');

figure(2)
%hold off;
axis( [0 400 200 1400]);
xlabel('number of X molecules');
ylabel('number of Y molecules');
title('Chemical system with SNIPER bifurcation');

%  tc{1}
%  tc{2}
%  tc{3}

%
%  plot( tc{1},Ac{1}(1,:), tc{2},Ac{2}(1,:), tc{3},Ac{3}(1,:), tc{4},Ac{4}(1,:), tc{5},Ac{5}(1,:), ...
%        tc{6},Ac{6}(1,:), tc{7},Ac{7}(1,:), tc{8},Ac{8}(1,:), tc{9},Ac{9}(1,:), tc{10},Ac{10}(1,:), ...
%        tm,Mm, 'k--');
%  hold on;
%  plot( tm,Mm, 'k--', 'LineWidth',3);
%  hold off;


%plot( tpl,Apl(1,:), t,A(1,:));
%axis 'auto x'

