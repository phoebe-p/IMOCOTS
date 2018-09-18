function tr = regular_pyramids_i(varargin)
% only want 3 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:regular_pyramids_u:TooManyInputs', ...
        'requires at most 1 optional input');
end

optargs = {54.7};
optargs(1:numvarargs) = varargin;
[char_angle] = optargs{:};

Lx = 1;
Ly = 1;
h = Lx*tand(char_angle)/2;

% make pyramids bigger

x = [0 Lx/2 Lx 0 Lx];
y = [0 Ly/2 0 Ly Ly];
z = [0 -h 0 0 0];
tri = delaunay(x,y);
tr = triangulation(tri,x(:),y(:),z(:));
