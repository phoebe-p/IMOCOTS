function tr = random_pyramids_i(varargin)
% only want 3 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 3
    error('myfuns:random_pyramids_u:TooManyInputs', ...
        'requires at most 3 optional inputs');
end
load('AFM.mat', 'pyramids_20um');
z = pyramids_20um;

optargs = {20 20 max(z(:))};
optargs(1:numvarargs) = varargin;
[x_width, y_width, z_max] = optargs{:};


z = -z_max*z/max(z(:));
[n_points_x, n_points_y] = size(z);
[x,y] = meshgrid(linspace(0, x_width, n_points_x),linspace(0, y_width, n_points_y));
tri = delaunay(x,y);
tr = triangulation(tri,x(:),y(:),z(:));
p = trisurf(tr);
p1 = reducepatch(p, 300);
pp = patch('Faces',p1.faces,'Vertices',p1.vertices);
tr = triangulation(pp.Faces,pp.Vertices);
