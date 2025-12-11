f1 = 25; f2 = 35; sx = 5 * 1.0e-3; sy = 5 * 1.0e-3; Cx = 1280/2; Cy = 1024/2;

K1 = [f1/sx, 0, Cx; 0, f1/sy, Cy; 0 0 1]; 
K2 = [f2/sx, 0, Cx; 0, f2/sy, Cy; 0 0 1]; 

R = [0 0 -1; 0 1 0; 1 0 0];
T = [500 * cosd(45); 0; 500 * sind(45)];

E = R' * myskew(T);
E1 = myskew(R' * T) * R';

F = inv(K2') * E * inv(K1);

det(E)
det(F)
svd(E)

e1 = K1 * T;
e2 = - K2 * R' * T;

M = [10; 20; 350];
m1 = K1 * M;
m1 = m1 / m1(3);
M_2 = inv([R, T; 0 0 0 1]) * [M; 1];
m2 = K2 * M_2(1:3);
m2 = m2 / m2(3);

m1_norm = inv(K1) * m1;
m2_norm = inv(K2) * m2;

ep2 = F * m1;
ep1 = F' * m2;

% Ppties
disp(m2' * F * m1)
disp(m2_norm' * E * m1_norm)

disp(F * e1)
disp(F' * e2)

disp(ep1' * m1)
disp(ep2' * m2)

disp(ep1 / norm(ep1) - cross(m1, e1) / norm(cross(m1, e1)))
disp(ep2 / norm(ep2) - cross(m2, e2) / norm(cross(m2, e2)))

F / F(3,3)
disp('End');
%% 
function res = myskew(v)
  res = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end