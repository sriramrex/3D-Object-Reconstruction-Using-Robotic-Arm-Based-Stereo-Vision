clear; clc; close all;

%% Read the images and assign them to the respective variables

folder1 = 'img/';
file1 = 'img';
extension1 = '.jpg';
filepath1 = [folder1,file1,extension1];
OriginalImg = imread(filepath1);
figure,imshow(OriginalImg)

%Image to gray

GrayImg = rgb2gray(OriginalImg);
figure, imshow(GrayImg);
plot(GrayImg(:,1000));

%% Gaussian Curve fitting

x = 1:1200;
x0 = 352;
sigma = 2;
y = 250*exp((-((x - x0) / sigma).^2) + 5 + 1 * randn(1, 1200));
figure
plot(x,y);
% First look for the maximum of the peak
[~, xMax] = max(y);
%%
mask = [-1 -2 -3 0 3 2 1];
yConv = conv(y, mask, 'same');
figure
plot(x, yConv);
intervalZCross = [xMax - 1: xMax + 2];
A = [x(intervalZCross)', yConv(intervalZCross)',ones(size(x(intervalZCross)'))];
[~, ~, V] = svd(A);
Recta = V(:, end);
hold on
plot(x(intervalZCross), (- Recta(1) * x(intervalZCross) - Recta(3)) /Recta(2), 'r');
x0ZCross = - Recta(3) / Recta(1);







