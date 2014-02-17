x = rand(4,1);
fprintf('vertices(1,:) = [%12.10f, %12.10f, %12.10f, %12.10f]\n', x)
y = rand(4,1);
fprintf('vertices(2,:) = [%12.10f, %12.10f, %12.10f, %12.10f]\n', y)
fprintf('vertices(3,:) = [0, 0, 0, 0]\n')
fprintf('area = %12.10f\n', polyarea(x,y))
vecn = rand(1,3);
vecn = vecn / sqrt(sum(vecn.*vecn));
fprintf('N = [%12.10f, %12.10f, %12.10f]\n', vecn(:))
alpha = rand(1,1)*360;
fprintf('alpha = %12.10f\n', alpha)

Rx = zeros(4,3);
for ipoint = 1:4
   r = [x(ipoint), y(ipoint), 0];
   Rx(ipoint,:) = vecn*dot(vecn,r) + cosd(alpha) * cross(cross(vecn,r),vecn) + sind(alpha) * cross(vecn, r);
end

fprintf('vertices_rot(1,:) = [%12.10f, %12.10f, %12.10f, %12.10f]\n', Rx(:,1))
fprintf('vertices_rot(2,:) = [%12.10f, %12.10f, %12.10f, %12.10f]\n', Rx(:,2))
fprintf('vertices_rot(3,:) = [%12.10f, %12.10f, %12.10f, %12.10f]\n', Rx(:,3))
