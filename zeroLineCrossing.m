clear all; close all; clc;
%%
x = 1:1200;
x0 = 512.7;
sigma = 10;
y = 250 * exp(-((x - x0) / sigma).^2) + 5 + 1 * randn(1, 1200);
[~, xMax] = max(y);
mask = [-1 -2 -3 0 3 2 1];
yConv = conv(y, mask, 'same');
figure
plot(x, yConv);

% We look at the zero crossing of this curve.
% First fitting a line a * x + b * y + c = 0
intervalZCross = [xMax - 4: xMax + 4];
A = [x(intervalZCross)', yConv(intervalZCross)', ones(size(x(intervalZCross)'))];
[~, ~, V] = svd(A);
Recta = V(:, end);
hold on
plot(x(intervalZCross), (- Recta(1) * x(intervalZCross) - Recta(3)) / Recta(2), 'r');
x0ZCross = - Recta(3) / Recta(1);