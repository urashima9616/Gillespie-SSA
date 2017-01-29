% degradation:  A --> 0 with rate constant k

clear all;

k = 0.1; % sec^(-1)

q = 1;
Nchs = 1;

knu = [ k ];
c = [1];
p = [0];
A0 = [20];
tfin = 30;
%    knu ... 1-by-q row vector; rate constants of the chemical reactions;
%            knu = k/nu^(order-1), where nu=volume and order=order of the chemical reaction
%    c ... Nchs-by-q; c(j,i) = coeficient at the j-th chemical species at the reactant side of the i-th reaction
%    p ... Nchs-by-q; p(j,i) = coeficient at the j-th chemical species at the product side of the i-th reaction
%    A0 ... Nchs-by-1 column vector;
%           A0(j) = number of molecules of j-th chemical species at the initial time t=0
%    tfin ... final time; if the simulation raches tfin it stops

% plot 10 realizations
figure(1)
clf;
fontsz = 20;
set(gca,'FontSize',fontsz);
%hold on;
for n = 1 : 10
  [A,t] = gillespiessa(knu, c,p, A0,tfin);
  [Apl,tpl] = mkplotdata(A,t);
  Ac{n} = Apl;
  tc{n} = tpl;
  %plot( tpl,Apl(1,:));
end % for n
%hold off;

%  tc{1}
%  tc{2}
%  tc{3}

tm = [0:0.1:30];
Mm = A0(1)*exp(-k.*tm);

plot( tc{1},Ac{1}(1,:), tc{2},Ac{2}(1,:), tc{3},Ac{3}(1,:), tc{4},Ac{4}(1,:), tc{5},Ac{5}(1,:), ...
      tc{6},Ac{6}(1,:), tc{7},Ac{7}(1,:), tc{8},Ac{8}(1,:), tc{9},Ac{9}(1,:), tc{10},Ac{10}(1,:), ...
      tm,Mm, 'k--');
hold on;
plot( tm,Mm, 'k--', 'LineWidth',3);
hold off;
axis( [0 tfin 0 A0(1)+1]);
text(18,18,['k=' num2str(k)], 'FontSize', fontsz);
xlabel('time');
ylabel('number of molecules');
title('Degradation A --> 0');

%plot( tpl,Apl(1,:), t,A(1,:));
%axis 'auto x'
