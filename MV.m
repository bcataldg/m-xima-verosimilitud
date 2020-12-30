function [mu sigma] = MV(VAR)

%var = [1xn]

m=length(VAR);
GM=1:m;

var=log(VAR);

p_O=@(O,mu,SIG) mvnpdf(O,mu,SIG);

med = nanmean(var);

des = nanstd(var);

myfun=@(O) -sum(log(p_O(var',O(1),O(2)^2)));

O0=[med des];   % semilla

A = []; b = []; % restricciones

[alpha,val,~,~,~,grad,~]=fmincon(myfun,O0,A,b);
mu=alpha(1);
sigma=alpha(2);

Pexc=1-normcdf(log(0.018),mu ,sigma); 

