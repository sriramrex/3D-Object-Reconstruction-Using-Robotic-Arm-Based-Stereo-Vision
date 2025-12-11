clear all; close all; clc;
%%
% Calibrate the laser camera with the cone
img = imread('img.jpg');
imshow(img);
load('H_matrix.mat');
invH = pinv(H_matrix);

% Black out the non interesting region
imgInteres = img * 0;

reg{1} = [600, 741];
reg{2} = [750, 780];
reg{3} = [788, 959];
reg{4} = [970, 995];
reg{5} = [1009, 1200];
reg{6} = [1205, 1230];
reg{7} = [1235, 1445];
reg{8} = [1456, 1473];
reg{9} = [1478, 1568];
reg{10} = [1581, 1627];
imgGray = rgb2gray(img);


for i = 1:length(reg)
    imgInteres(:, reg{i}(1):reg{i}(2), :) = ...
        img(:, reg{i}(1):reg{i}(2), :);
end

% figure;
% imshow(imgInteres);
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

% figure;
% imshow(imgMax);

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

for i = 1:length(reg)-1
    PointsIm(:, i) = -pinv([Recta{i}(1:2);Recta{i + 1}(1:2)])*[Recta{i}(3); Recta{i + 1}(3)];
    %hold on
    plot(PointsIm(1, i), PointsIm(2, i), 'ob');
end
hold off





%%
%Add randomness to the points of cone
% randPtsCone = [];
% sigmaPt = 0.1;
% for m = 1:101
%     randPtsCone= [PtsConeI(1,m)+sigmaPt.*randn(100,1),PtsConeI(2,m)+sigmaPt.*randn(100,1)];
%     plot(randPtsCone(:,1),randPtsCone(:,2),'ok');
% end

%%

PtsCono = [75 75 62 62 49 49 36 36 23 23 -23; ...
    0 50 50 100 100 150 150 200 200 222 222];
PtsConeCorner = H_matrix * [PtsCono; ones(1, size(PtsCono, 2))];
PtsConeCorner1 = PtsConeCorner ./ PtsConeCorner(3, :);
% figure;
% plot(PtsConeCorner1(1:2,1),PtsConeCorner1(1:2,2),'or');
% hold on

for region = 1:length(reg)
    ptsRegI = [];
    LineReg = cross(PtsConeCorner1(:,region), PtsConeCorner1(:, region + 1));
    for i = reg{region}(1) : reg{region}(2)
        ptsI = -(LineReg(1) * i + LineReg(3)) / LineReg(2);
        ptsRegI = [ptsRegI, [i;ptsI]];
    end

    ptsConereg = invH * [ptsRegI; ones(1, size(ptsRegI, 2))];
    ptsConereg = ptsConereg ./ ptsConereg(3,:);

    ptsConeRegs{region} = ptsConereg;
end

figure()
for i = 1:length(reg)
    ptsCornRegI = ptsConeRegs{i};
    plot(ptsCornRegI(1,:), ptsCornRegI(2,:), '.'); axis equal;
    hold on
end
 
%% Randomness

% Add randomness to the give points 
randPtsConeRight = [];
randPtsConeLeft = [];
sigmaPt = 0.01;
Error = [];

for n = 1:length(reg)
   randPtsConeRight = [];
   randPtsConeLeft = [];
    if mod(n,2)==1
        p = length(ptsConeRegs{n});
        for i1 =  1:p
            current_pts = ptsConeRegs{n};
            randPt = sigmaPt.*randn(1,1);
            randConeRight= [current_pts(1,i1)+randPt,current_pts(2,i1)];
            randConeLeft= [current_pts(1,i1)-randPt,current_pts(2,i1)];
            plot(randConeRight(:,1),randConeRight(:,2),'ok');
            %plot(randConeLeft(:,1),randConeLeft(:,2),'ok'); 
            Error= [Error,sqrt((randConeRight(:,1)-randConeLeft(:,1))^2+(randConeRight(:,2)-randConeLeft(:,2))^2)];
            randPtsConeRight = [randPtsConeRight;randConeRight];
            randPtsConeLeft = [randPtsConeLeft;randConeLeft];
        end
      

    else
        p = length(ptsConeRegs{n});
        for i2 =  1:p
            current_pts = ptsConeRegs{n};
             randPt = sigmaPt.*randn(1,1);
            randConeRight= [current_pts(1,i2),current_pts(2,i2)+randPt];
            randConeLeft= [current_pts(1,i2),current_pts(2,i2)-randPt];
             plot(randConeRight(:,1),randConeRight(:,2),'ok');
            %plot(randConeLeft(:,1),randConeLeft(:,2),'ok');
            Error= [Error,sqrt((randConeRight(:,1)-randConeLeft(:,1))^2+(randConeRight(:,2)-randConeLeft(:,2))^2)];
            randPtsConeRight = [randPtsConeRight;randConeRight];
            randPtsConeLeft = [randPtsConeLeft;randConeLeft];
        end
      
    
    end
   RandomOne{n}=[randPtsConeRight'];
   RandomTwo{n}=[randPtsConeLeft'];

end

%% Reproject the points on the image plane


figure()

for i = 1:length(reg)
    ptsCornRegI = ptsConeRegs{i};
    RandOne = RandomOne{i};
    RandTwo = RandomTwo{i};
    Image_newPts = H_matrix*ptsCornRegI;
    RandomOne_Img = H_matrix*[RandOne;ones(1,length(RandOne))];
    RandomTwo_Img = H_matrix*[RandTwo;ones(1,length(RandTwo))];
    plot(Image_newPts(1,:), Image_newPts(2,:), '.'); axis equal;
   
    plot(RandomOne_Img(1,:), RandomOne_Img(2,:), 'or'); axis equal;
   
    %plot(RandomTwo_Img(1,:), RandomTwo_Img(2,:), 'ob'); axis equal;
    hold on
    ImageRand{i} = [RandomOne_Img];
end

%% Randomness_Image

% Add randomness to the give  image points 

sigmaImg = 0.001;

for n = 1:length(reg)
    randImgCone = [];
    q = length(ImageRand{n});
    for i3 =  1:q
        RandImgPts = ImageRand{n};
        RandImgOff = sigmaImg.*randn(1,1);
        RadPts = [ RandImgPts(1,i3), RandImgPts(2,i3)+RandImgOff];

        randImgCone = [randImgCone; RadPts ];
    end
    RandomImg{n}=[randImgCone'];
end




for i = 1:length(reg)
    ptsCornRegI = ptsConeRegs{i};
    RandImg = RandomImg{i};
    plot(RandImg(1,:), RandImg(2,:), 'ob'); axis equal;
    hold on
end

%% combing the points set in the image:
ErrorDist = [];
for i = 1:length(reg)
     q = length(ImageRand{i});
     Holder = [];
     dist = [];
     Rand1 = ImageRand{i}(1:2,:);
     Rand2 = RandomImg{i};
    for j1 = 1:q
    Holder = [Holder,Rand1(:,j1)];
    Holder = [Holder,Rand2(:,j1)];
    dist = [dist,abs(Rand1(2,j1)-Rand2(2,j1))];
    end
    PtsSet{i} = [Holder];
    ErrorDist{i} = [dist];
end









%%
% Perform the fitting of the lines
for i = 1:length(reg)
  rang = PtsSet{i}(1,:);
  Sub = PtsSet{i}(2,:);
  Recta{i} = fitLine(rang, Sub);
end
hold on 
% Plot the lines on the image
for i = 1:length(reg)
    lt = min(PtsSet{i}(1,:))-0.01;
    rt = max(PtsSet{i}(1,:))+0.01;
  intervalX = [lt,PtsSet{i}(1,:),rt];
  y = -(Recta{i}(1) * intervalX + Recta{i}(3)) / Recta{i}(2);
  plot(intervalX, y, '-g');
end

%%
% Calculate the intersection of the lines
for i = 1:length(reg) - 1
 PointsIm(:, i) = -pinv([Recta{i}(1:2);Recta{i + 1}(1:2)])*[Recta{i}(3); Recta{i + 1}(3)];
 hold on
 plot(PointsIm(1, i), PointsIm(2, i), 'ok');
end

%%
%hold off
% Calculate the homograhy     
PtsCono = [75 62 62 49 49 36 36 23 23; ...
           50 50 100 100 150 150 200 200 222];

ErrorCone = [];
     
H_matrixNew = get_homography(PointsIm, PtsCono);
NewInvH = pinv(H_matrixNew);


for re= 1:length(reg)
    ptsImgreg = NewInvH * [ImageRand{re}(1:2,:); ones(1, length(ImageRand{re}))];
    ptsImgreg = ptsImgreg ./ ptsImgreg(3,:);
    PtsImgRegs{re} = [ptsImgreg];
end


hold on
for i = 1:length(reg)
    ptsImgRegI = PtsImgRegs{i};
    ptsCornRegI = RandomOne{i};
    plot(ptsImgRegI(1,:), ptsImgRegI(2,:), '.b'); axis equal;
    hold on
    ErrorCone{i}= [sqrt((ptsCornRegI(1,:)-ptsImgRegI(1,:)).^2+(ptsCornRegI(2,:)-ptsImgRegI(2,:)).^2)];
end

%%
%Plot Mean and Standard Deviationn

x1 = ErrorCone{1}; % Centered and sigma = 1
ystd = std(x1);
ymean = mean(x1);




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
