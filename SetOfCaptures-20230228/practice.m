imgInteres = img *0;

reg{1} = [695, 1025]; 

imgGray = rgb2gray(img);

for i = 695:1025
  imgInteres(:, i, :) = ...
     img(:, i, :);
end

figure;
imshow(imgInteres);

%% 
% Calculate the subpixel peaks
for i = 695:1025
  %Calculate the max peaks
 [~, indexes{i}] = max(imgGray(:, i, 1));

  % Calcular los subpeaks
   delta = 2;
   rang = i;
 initMax = indexes{i};
   indexSub{i} = maxSubPeak(imgGray, rang, initMax, delta);
end

imgMax = img;
ptPlot = [];
colorValues =[];
for i = 695:1025
  %for j = 1:length(indexes{i})
    imgMax(indexes{i},i, 1) = 0;
    imgMax(indexes{i},i, 2) = 255;
    imgMax(indexes{i},i, 3) = 0;

    indexesSub = abs(round(indexSub{i}));
    imgMax(indexesSub,i, 1) = 255;
    imgMax(indexesSub,i, 2) = 255;
    imgMax(indexesSub, i, 3) = 0;
    ptPlot(i-694,:) =[indexesSub,i] ;
    
  %end
end
x_points = ptPlot(:,1);
y_points = ptPlot(:,2);

for i = 1:length(x_points)
colorValues(i,:) = [R(y_points(i),x_points(i)),G(y_points(i),x_points(i)),B(y_points(i),x_points(i))];
end


figure;
imshow(imgMax);
hold on
plot(ptPlot(:,2),ptPlot(:,1),'ob')
%[fitresult, gof] = curveFit(y_points, x_points);

%%
PtsTourn = [];
for k = 1:331
 ptsConeI = pinv(H_matrix) * [ptPlot(k,2); ptPlot(k,1); 1];
  ptsConeI = ptsConeI / ptsConeI(3);
  
  if (ptsConeI(2) < 0) % Points with Z positive
    continue;
  end
  
  PtsTourn = [PtsTourn, ptsConeI(1:2, :)];
  %RGBI = imgNolaser(round(indexSubI), i, :);
  %ColorTourn = [ColorTourn, [RGBI(:, :,1); RGBI(:, :,2); RGBI(:, :,3)]];
end
% i = 1;
% 
% obj = zeros(3,1065);
% for z=1:360/1065:360
% rz = rotz(z);
% obj(:,i) = rz*[PtsTourn(1,i);PtsTourn(2,i);1];
% i = i+1;
% end

figure;
 plot(PtsTourn(1, :), PtsTourn(2, :), 'ok'); axis equal


%%
% Calculate the cone points with the total detected points
PtsTourn = [];
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
end

i = 1;
obj = zeros(3,length(PtsTourn));
for z=1:360/length(PtsTourn):360
rz = rotz(z);
obj(:,i) = rz*[PtsTourn(1,i);PtsTourn(2,i);1];
i = i+1;
end

figure
plot(PtsTourn(1, :), PtsTourn(2, :), 'ok'); axis equal


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

