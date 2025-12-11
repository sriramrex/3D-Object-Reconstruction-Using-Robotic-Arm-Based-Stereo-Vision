clear all; close all; clc;
%%
% Calibrate the laser camera with the cone
load('cameraParams');
load('H_matrix.mat');
directory = 'Captures_Cup';
Files = dir(strcat(directory, '\imgLaser*.jpg'));
FilesColor = dir(strcat(directory, '\imgNoLaser*.jpg'));
file = fopen(strcat(directory, '\cupPts.txt'), 'w');
for i0 =1:length(Files)
    disp(strcat(directory, '\', Files(i0).name));
    img = imread(strcat(directory, '\', Files(i0).name));
    img = undistortImage(img, cameraParams);
    imgColor = imread(strcat(directory, '\', FilesColor(i0).name));
    imgColor = undistortImage(imgColor, cameraParams);
    %imshow(img);

    imgGray = rgb2gray(img);
    PtsTourn = [];
    PtsColor = [];

    R = imgColor(:,:,1);
    G = imgColor(:,:,2);
    B = imgColor(:,:,3);
    %end
    %%
    % Black out the non interesting region
    imgInteres = img *0;

    reg{1} = [1, size(imgGray, 2)];

    imgGray = rgb2gray(img);

    for i = 1:size(imgGray, 2)
        imgInteres(:, i, :) = ...
            img(:, i, :);
    end

    %figure;
    %imshow(imgInteres);

    % Calculate the subpixel peaks
    for i = 1:size(imgGray, 2)
        %Calculate the max peaks
        [~, indexes{i}] = max(imgGray(:, i, 1));

        % Calcular los subpeaks
        delta = 2;
        rang = i;
        initMax = indexes{i};
        %indexSub = [];
        indexSub{i} = maxSubPeak(imgGray, rang, initMax, delta);
    end

    imgMax = img;
    ptPlot = zeros(size(imgGray, 2), 2);

    colorValues =[];
    for i = 1:size(imgGray, 2)

        %for j = 1:length(indexes{i})
        if (isempty(indexes{i}) || indexes{i} < 1 || ...
                indexes{i} > size(imgGray, 1))
            continue;
        end

        indexesSub = round(indexSub{i});

        if (indexesSub < 1 || indexesSub > size(imgGray, 1))
            continue;
        end

        imgMax(indexes{i},i, 1) = 0;
        imgMax(indexes{i},i, 2) = 255;
        imgMax(indexes{i},i, 3) = 0;

        imgMax(indexesSub,i, 1) = 255;
        imgMax(indexesSub,i, 2) = 255;
        imgMax(indexesSub,i, 3) = 0;
        ptPlot(i,:) =[indexesSub,i] ;

    end
    x_points = ptPlot(:,1);
    y_points = ptPlot(:,2);



    PtsTourn = [];
    colorValues = [];

    for k = 1:size(imgGray, 2)
        if (ptPlot(k, 2) == 0 || ptPlot(k, 1) == 0)
            continue;
        end

        ptsConeI = pinv(H_matrix) * [ptPlot(k,2); ptPlot(k,1); 1];
        ptsConeI = ptsConeI / ptsConeI(3);

        if (ptsConeI(2) < 0 || ptsConeI(1) < 0) % Points with Z positive
            continue; 
        end

        PtsTourn = [PtsTourn, ptsConeI(1:2, :)];


        x = round(x_points(k, 1));
        y = round(y_points(k, 1));
        colorValues = [colorValues; [R(x, y),G(x,...
                          y),B(x, y)]];

    end

    theta = i0 * 1.8; % in degrees
    for j = 1:size(PtsTourn, 2)
        PtsIJ = PtsTourn(1:2, j);
        PtsXYZIJ = [PtsIJ(1) * cosd(theta); PtsIJ(1) * sind(theta); PtsIJ(2)];
        fprintf(file, '%f, %f, %f, %f, %f, %f\n', PtsXYZIJ(1), PtsXYZIJ(2), PtsXYZIJ(3),...
            colorValues(j,1), colorValues(j,2), colorValues(j,3));
    end
end
fclose(file);

%%
function [peak] = maxSubPeak(imgGray, rang,initMax,delta)
% rang = 894;
% initMax = indexes{894};
imgt = imgGray(:, rang);
% figure,imshow(imgt);
% figure;
for j = 1:length(rang)
    %plot(imgt(:,j));
    line = imgt(:,j);
    x0 = initMax(j);

    if (x0 < delta || x0 > size(imgGray, 1) - delta)
        peak(j) = 0;
        continue;
    end

    if (line(x0) < 100)
        peak(j) = 0;
        continue;
    end

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
