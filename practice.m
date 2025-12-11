%%


% PointsIm1 = H_matrix * [PointsIm; ones(1, size(PointsIm, 2))];
% PointsIm1 = PointsIm1 ./ PointsIm1(3, :);
% figure;
% plot(PointsIm1(1:2,1),PointsIm1(1:2,2),'-r');
% hold on
% c = 1;
% PtsConeI = zeros(3,101);
% for n = 8625.62:43.8598:13011.6
%     x1 = n;
%     y1 = 5000;
%     x2 = n;
%     y2 = 9000;
%     yi = (PointsIm1(1,1)-PointsIm1(2,1))*(y1-y2)-(PointsIm1(1,2)-PointsIm1(2,2))*(x1-x2);
%     xx = (PointsIm1(1,1)*PointsIm1(2,2)-PointsIm1(1,2)*PointsIm1(2,1))*(x1-x2)-(PointsIm1(1,1)-PointsIm1(2,1))*(x1*y2-y1*x2);
%     xy = (PointsIm1(1,1)*PointsIm1(2,2)-PointsIm1(1,2)*PointsIm1(2,1))*(y1-y2)-(PointsIm1(1,2)-PointsIm1(2,2))*(x1*y2-y1*x2);
%     pxi = xx/yi;
%     pyi = xy/yi;
%     plot([x1,x2],[y1,y2],'-g');
%     hold on
%     plot(pxi,pyi,'ob');
%     % Calculate the cone points with the inverse of the homography
%     PtsCone1 = invH * [pxi;pyi; 1];
%     PtsConeI(:,c) = [PtsCone1 ./ PtsCone1(3, :)];
%     c=c+1;
% end
% 
% 
% hold off

%%
%Projected points using homography inverse
% figure;
% plot(PtsConeI(1,:),PtsConeI(2,:),'or');
% 
% hold on;

%plot(ptsConeCorners(1,:),ptsConeCorners(2,:),'.k'); axis equal