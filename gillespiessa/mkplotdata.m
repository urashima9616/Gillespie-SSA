function [Apl,tpl] = mkplotdata(A,t)
% Prepare data for stairs-like graphs.
% INPUT
% A,t ... output data from gillespiessa.m
% OUTPUT
% Apl,tpl ... a modification of A and t. If Apl, tpl is used with standard 
%     Matlab plot commands then a stairs-like graphs are produced. Such 
%     graphs represent functions which are constant over each time step,
%     which is the correct interpretation of the output of the Gillespie SSA.


Nchs = size(A,1);

N = length(t);
Npl = 2*(N-2) + 2;

Apl = zeros(Nchs,Npl);
for m = 1:N-1
 Apl(:,2*m-1) = A(:,m);
 Apl(:,2*m)   = A(:,m);
end % for

tpl = zeros(1,Npl);
tpl(1) = t(1);
tpl(Npl) = t(N);
for m = 2:N-1
 tpl(2*m-2) = t(m);
 tpl(2*m-1) = t(m);
end % for
