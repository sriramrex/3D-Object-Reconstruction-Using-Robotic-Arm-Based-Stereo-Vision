function [peak] = maxSubPeak(imgGray, rang,initMax,delta)

imgt = imgGray(:, rang);
% figure,imshow(imgt);
% figure;
for j = 1:length(rang)
plot(imgt(:,j));
line = imgt(:,j);
x0 = initMax(j);

interval = [x0 - delta: x0 + delta];
logy = log(abs(double(line(interval))) + 1)'; % to prevent log(0)
xi = double(interval)';
B = logy;
A = [-xi.^2, 2 * xi, - ones(size(xi))];
Res = pinv(A) * B';
sigma2Res = 1 / Res(1);
x0Exp = Res(2) * sigma2Res;


peak(j) = x0Exp;
end
end