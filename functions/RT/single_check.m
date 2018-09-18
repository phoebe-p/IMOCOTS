function[results] = single_check(tr, FN, r_a, d, side)
results = [];
triangle_coords = tr.Points(tr.ConnectivityList',:); % x y and z coordinates for corners of triangles.
% every group of 3 rows are three points corresponding to one triangle on 
% the texture, the columns are x y and z coords.

P_0s = triangle_coords(1:3:end,:)'; % each column is one the corners of the triangles (3D coordinate)
P_1s = triangle_coords(2:3:end,:)'; % each column is the second corner
P_2s = triangle_coords(3:3:end,:)'; % each column is the third corner

% below I am calculating t, u, v. 
% see https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
% Point of intersection is at: r_a + t*d
% general point on the triangle: P0 + (P1-P0)*u + (P2-P0)*v
% where P0, P1 and P2 are the 3D coordinates of the corner of a specific
% triangle.
% t has to be > 0, otherwise the intersection occurs 'behind' the ray.
% u + v has to be <= 1, otherwise the intersection occurs in the plane the 
% triangle is in but outside the bounds of the triangles
[~, width] = size(P_0s);
crossP = cross(P_1s-P_0s, P_2s - P_0s); 
pref = 1./dot(repmat(-d, 1, width), crossP, 1);
corner = r_a-P_0s;
t = pref.*dot(crossP, corner, 1);
u = pref.*dot(cross(P_2s-P_0s, repmat(-d, 1, width)), corner, 1);
v = pref.*dot(cross(repmat(-d, 1, width), P_1s-P_0s), corner, 1);
As = [t; u; v];

which_intersect = (As(2,:) + As(3,:)) <= 1 & all(As(2:3,:) >= 0) & As(1,:) > 0;
t(which_intersect == 0) = NaN;
[t, ind] = min(t);
%disp('Intersection!')
intersn = r_a + t*d;
N = side*FN(ind,:)';
ThetaInDegrees = atan2d(norm(cross(N,-d)),dot(N,-d));
results.intersn = intersn;
results.theta = ThetaInDegrees;
results.N = N;

% result.intersn will be NaNs if no intersection

end

