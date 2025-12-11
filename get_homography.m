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
A(i,:) = [z_t*x, z_t*y, z_t*z, 0, 0, 0, -x_t*x, -x_t*y, -x_t*z];
A(j,:) = [0, 0, 0, z_t*x, z_t*y, z_t*z, -y_t*x, -y_t*y, -y_t*z];
end

[U,S,V] = svd(A);
 
H_matrix  = V(:,end);
H_matrix = reshape(H_matrix,3,3)';

end