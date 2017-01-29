% production and degradation:
%   A --> 0 with rate constant k1
%   0 --> A with rate constant k2

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

k1 = 0.1; % sec^(-1)
k2nu = 1; % sec^(-1)

q = 2;
Nchs = 1;

knu = [ k1 k2nu ];
c = [1, 0];
p = [0, 1];
A0 = [0];
tfin = 100;

%ode45

% plot Nreal realizations
Nreal = 5;
figure(1)
clf;
fontsz = 20;
set(gca,'FontSize',fontsz);
hold on;
for n = 1 : Nreal
  [A,t] = gillespiessa(knu, c,p, A0,tfin);
  [Apl,tpl] = mkplotdata(A,t);
  Ac{n} = Apl;
  tc{n} = tpl;
  co = get(gca,'ColorOrder');
  plot( tpl,Apl(1,:), 'Color', co(mod(n,7)+1,:));
end % for n

tm = [0:0.1:tfin];
Mm = k2nu/k1*(1-exp(-k1.*tm));
plot( tm,Mm, 'k--', 'LineWidth',3);

hold off;

text(3,19,['k1=' num2str(k1)], 'FontSize', fontsz);
text(3,17,['k2=' num2str(k2nu)], 'FontSize', fontsz);
xlabel('time');
ylabel('number of molecules');
title('Production and degradation A <--> 0');

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
axis( [0 tfin 0 21]);

%plot( tpl,Apl(1,:), t,A(1,:));
%axis 'auto x'
