clear all; close all; clc;
% Calibrate the laser camera with the cone

% Take image of cone
% cam = webcam(1);
% img = snapshot(cam);

% For debugging
%imwrite(img, 'img.jpg');
img = imread('img.jpg');
imshow(img);

% Black out the non interesting region
imgInteres = img * 0;

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
   indexSub{i} = maxSubPeak(imgGray, rang, initMax, delta);
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
% Perform the fitting of the lines
for i = 1:length(reg)
  rang = reg{i}(1):reg{i}(2);
  Sub = indexSub{i};
  Recta{i} = fitLine(rang, Sub);
end

hold on 
% Plot the lines on the image
for i = 1:length(reg)
  intervalX = reg{i}(1, 1) - 20:reg{i}(1, 2) + 20;
  y = -(Recta{i}(1) * intervalX + Recta{i}(3)) / Recta{i}(2);
  plot(intervalX, y, '-b');
end

%%
% Calculate the intersection of the lines
% for i = 1:length(reg) - 1
%   PointsIm(:, i) = -pinv([Recta{i}(1:2)'; Recta{i + 1}(1:2)']) * ...
%     [Recta{i}(3); Recta{i + 1}(3)];
%   hold on
%   plot(PointsIm(1, i), PointsIm(2, i), 'oy');
% end

for i = 1:length(reg) - 1
 PointsIm(:, i) = -pinv([Recta{i}(1:2);Recta{i + 1}(1:2)])*[Recta{i}(3); Recta{i + 1}(3)];
 hold on
 plot(PointsIm(1, i), PointsIm(2, i), 'ob');
end
%%
% Calculate the homograhy     
PtsCono = [75 62 62 49 49 36 36 23 23; ...
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

% Calculate the cone points with the total detected points
PtsTourn = [];
%ColorTourn = [];
%imgNolaser = imread('imgNolaser.jpg');

%%
for i = 1:size(imgGray, 2)
  [peakMax, indexesI] = ...
    max(imgGray(:, i));
  
  if (peakMax < 150)
    continue;
  end
  
  delta = 7;
  indexSubI = maxSubPeak(imgGray, i, indexesI, delta);
  
  ptsConeI = pinv(H_matrix) * [i; indexSubI; 1];
  ptsConeI = ptsConeI / ptsConeI(3);
  
  if (ptsConeI(2) < 0) % Points with Z positive
    continue;
  end
  
  PtsTourn = [PtsTourn, ptsConeI(1:2, :)];
  %RGBI = imgNolaser(round(indexSubI), i, :);
  %ColorTourn = [ColorTourn, [RGBI(:, :,1); RGBI(:, :,2); RGBI(:, :,3)]];
end
i = 1;

obj = zeros(3,1065);
for z=1:360/1065:360
rz = rotz(z);
obj(:,i) = rz*[PtsTourn(1,i);PtsTourn(2,i);1];
i = i+1;
end

figure;
 plot(PtsTourn(1, :), PtsTourn(2, :), 'ok'); axis equal

%  X = PtsTourn(1, :);
%  Y = PtsTourn(2, :);
% [X,Y,Z] = cylinder();
% surf(X,Y,222)



% Tourn the points 1.8 degrees 200 times to make a tour
 file = fopen('conoPts.txt', 'w');

%%
for i = 0:199
  theta = i * 1.8; % in degrees
  for j = 1:size(PtsTourn, 2)
    PtsIJ = PtsTourn(1:2, j);
    PtsXYZIJ = [PtsIJ(1) * cosd(theta); PtsIJ(1) * sind(theta); PtsIJ(2)];
%     fprintf(file, '%f, %f, %f, %f, %f, %f\n', PtsXYZIJ(1), PtsXYZIJ(2), PtsXYZIJ(3), ...
%       ColorTourn(1, j), ColorTourn(2, j), ColorTourn(3, j));
   %fprintf(file, '%f, %f, %f, %f, %f, %f\n', PtsXYZIJ(1), PtsXYZIJ(2), PtsXYZIJ(3));
    fprintf(file, '%f, %f, %f\n', PtsXYZIJ(1), PtsXYZIJ(2), PtsXYZIJ(3));
  end
end
% 
 fclose(file);
%%
function [H_matrix] = get_homography(PointsIm, PtsCono)

a1 = zeros(18);
for i = 1:9
x = PtsCono(1,i);
y = PtsCono(2,i);
z = 1;
x_t = PointsIm(1,i);
y_t = PointsIm(2,i);
z_t = 1;
j = i*2;
A(j - 1,:) = [z_t*x, z_t*y, z_t*z, 0, 0, 0, -x_t*x, -x_t*y, -x_t*z];
A(j,:) = [0, 0, 0, z_t*x, z_t*y, z_t*z, -y_t*x, -y_t*y, -y_t*z];
end

[U,S,V] = svd(A);

H_matrix = V(:, end);
H_matrix  = reshape(H_matrix, 3, 3)';

end
%%
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
%%
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