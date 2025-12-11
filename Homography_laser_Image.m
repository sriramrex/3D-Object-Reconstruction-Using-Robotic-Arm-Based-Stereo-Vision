clear; clc; close all;

%% Read the images and assign them to the respective variables

folder1 = 'img/';
file1 = 'img';
extension1 = '.jpg';
filepath1 = [folder1,file1,extension1];
img = imread(filepath1);
figure,imshow(img)

%% Black out the non interesting region
imgInteres = img*0;
reg{1} = [600, 741];
reg{2} = [750, 780];
reg{3} = [788, 959];
reg{4} = [970, 995];
reg{5} = [1009, 1200];
reg{6} = [1205, 1232];
reg{7} = [1235, 1445];
reg{8} = [1456, 1473];
reg{9} = [1478, 1568];
reg{10} = [1581, 1627];
imgGray = rgb2gray(img);
for i = 1:length(reg)
 imgInteres(:, reg{i}(1):reg{i}(2), :) = ...
 img(:, reg{i}(1):reg{i}(2), :);
end
figure;
imshow(imgInteres);

%%
% Calculate the subpixel peaks
for i = 1:length(reg)
 %Calculate the max peaks
 [~, indexes{i}] = max(imgGray(:, reg{i}(1):reg{i}(2), 1));
 % Calcular los subpeaks
 delta = 2;
 rang = reg{i}(1):reg{i}(2);
 initMax = indexes{i};
 indexSub{i} = maxSubPeak(imgGray,rang,initMax,delta);
end

imgMax = img;
for i = 1:10
 for j = 1:length(indexes{i})
 imgMax(indexes{i}(j), reg{i}(1) + j - 1, 1) = 0;
 imgMax(indexes{i}(j), reg{i}(1) + j - 1, 2) = 255;
 imgMax(indexes{i}(j), reg{i}(1) + j - 1, 3) = 0;
 indexesSub = round(indexSub{i}(j));
 imgMax(indexesSub, reg{i}(1) + j - 1, 1) = 255;
 imgMax(indexesSub, reg{i}(1) + j - 1, 2) = 255;
 imgMax(indexesSub, reg{i}(1) + j - 1, 3) = 0;
 end
end
figure;
imshow(imgMax);
%%
%Perform the fitting of the lines
for i = 1:length(reg)
rang = reg{i}(1):reg{i}(2);
initSub = indexSub{i};
 Recta{i} = fitLine(rang, initSub);
end
hold on
%Plot the lines on the image
for i = 1:length(reg)
 intervalX = reg{i}(1, 1) - 20:reg{i}(1, 2) + 20;
 y = -(Recta{i}(1) * intervalX + Recta{i}(3)) / Recta{i}(2);
 plot(intervalX, y, '-b');
end

%%
for i = 1:length(reg) - 1
 PointsIm(:, i) = -pinv([Recta{i}(1:2);Recta{i + 1}(1:2)])*[Recta{i}(3); Recta{i + 1}(3)];
 hold on
 plot(PointsIm(1, i), PointsIm(2, i), 'ok');
end
%%
PtsCono = [75 62 62 49 49 36 36 23 23;...
    50 50 100 100 150 150 200 200 222];

H_matrix = get_homography(PointsIm, PtsCono);

% Verify the homography
PointsIm1 = H_matrix * [PtsCono; ones(1, size(PtsCono, 2))];
PointsIm1 = PointsIm1 ./ PointsIm1(3, :);

Err = PointsIm - PointsIm1(1:2, :);

% Calculate the cone points with the inverse of the homography
PtsCone1 = pinv(H_matrix) * [PointsIm; ones(1, size(PointsIm, 2))];
PtsCone1 = PtsCone1 ./ PtsCone1(3, :);

ErrCone = PtsCono - PtsCone1(1:2, :);

%%
for k = 1:9
x_t = PointsIm(1,k);
y_t = PointsIm(2,k);
z_t = 1;
Newpts(k,:) = [x_t,y_t,z_t]*H_matrix;
end
hold on
figure,plot(Newpts,'or');

