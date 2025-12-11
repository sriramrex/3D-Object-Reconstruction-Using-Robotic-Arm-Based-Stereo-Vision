function [coeff] = fitLine(rang,initSub )

for k = 1:length(rang)-1
y1 = initSub(k);
y2 = initSub(k+1);
x1 = rang(k);
x2 = rang(k+1);

a(k)= y1-y2;
b(k)= x2-x1;
c(k)= y2*x1-x2*y1;

end
coeff(1) = mean(a);
coeff(2) = mean(b);
coeff(3) = mean(c);

end