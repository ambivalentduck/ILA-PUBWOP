function analyze()

load schizo.csv

global A B C D

A=schizo(:,1:end-2);
C=schizo(:,end-1);
B=C>0;
B(~C)=1e10; %Infinite control cost for uncontrollable things
D=schizo(:,end);

C=C/10;
D=D/10;

f=find(B);

starting_eig=eig(A)
starting_cost=cost(zeros(size(B)))

b=fminunc(@cost,zeros(size(B)))

final_cost=cost(b)
e=eig(A+diag(b.*B))

function out=cost(in)
global A B C D

e=eig(A+diag(in.*B));
e(e<1)=0;

out=sum(D.*e.^2)+sum(C.*(in).^2);