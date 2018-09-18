function [Rmat, Tmat, Amat, angle_vector, theta_summary] = RT_Fresnel(incidence, behind, texture, options, wavelength)

% if ischar(texture)
%     [tr, options] = make_texture(texture, options);
% else
%     tr = texture;
% end
tr = texture;
usemean = options.usemean;

%figure
%trisurf(tr)

n_angle_bins = options.n_angle_bins;
theta_max = options.theta_max;
phi_max = options.phi_max;
c_azimuth = options.c_azimuth;
n_x = options.n_scans;
n_y = options.n_scans;
n_ray_angles = ceil(options.n_rays/(n_x*n_y));

% get optical constants
nk=load(incidence);
n1=interp1(nk(:,1)*1000,nk(:,2),wavelength)+1i*interp1(nk(:,1)*1000,nk(:,3),wavelength);
nk=load(behind);
n2=interp1(nk(:,1)*1000,nk(:,2),wavelength)+1i*interp1(nk(:,1)*1000,nk(:,3),wavelength);

[matrix_res, angle_vector, theta_summary] = ray_trace_3D(tr, usemean, n_ray_angles, n_angle_bins, n_x, n_y, theta_max, phi_max, c_azimuth, n1, n2);
% surface, n_rays, n_angle_bins, n_x,
% n_y, theta_max, phi_max, c_azimuth, n1, n2

% binning is out_loc, in_loc

Rmat = matrix_res(1:(length(angle_vector)/2), 1:(length(angle_vector)/2));
Tmat = flipud(matrix_res((length(angle_vector)/2 + 1):(end-1), 1:(length(angle_vector)/2)));
Amat = matrix_res(end, 1:length(angle_vector)/2);
end

function[matrix_res, angle_vector, theta_summary] = ray_trace_3D(tr, usemean, varargin)

numvarargs = length(varargin);
if numvarargs > 9
    error('ray_trace_3D takes at most 9 optional arguments');
end

% set defaults for optional inputs
optargs = {1000, 100, 10, 10, 180, 45, 0.25, 1, 3.14};
% n_rays, n_angle_bins, n_x, n_y, theta_max, phi_max, c_azimuth, n1, n2

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[n_rays, n_angle_bins, n_x, n_y, theta_max, phi_max, c_azimuth, n_1, n_2] = optargs{:};

FN = faceNormal(tr);
sin_a_b = linspace(0, 1, n_angle_bins+1); % number of bins is between 0 and 90
% will have the same number of
% bins between 90 and 180
% phi_max_geometry = 45;
theta_intervals=[asind(sin_a_b)  180-fliplr(asind(sin_a_b(1:end-1)))];
%theta_intervals = theta_intervals(theta_intervals <= theta_max);
theta_middle = mean([theta_intervals(1:end-1);theta_intervals(2:end)]);
phi_intervals = cell(length(theta_middle), 1);
angle_vector = [];
for j1 = 1:length(theta_middle)
    if theta_middle(j1) > 90
        ind = length(theta_intervals)-j1;
    else
        ind = j1;
    end
    phi_intervals{j1} = linspace(0, 90, ceil(c_azimuth*ind)+1);
    phi_middle = mean([phi_intervals{j1}(1:end-1);phi_intervals{j1}(2:end)]);
    angle_vector = [angle_vector; repmat(j1, length(phi_middle), 1) ...
        repmat(theta_middle(j1), length(phi_middle), 1) phi_middle'];
end

% can use angle_vector(:,1) to find where the theta entry start, then check
% the phi index. Then index is theta_index + phi_index - 1

% for the results, the length is the total number of rays checked, the
% number of NaNs is the absorption and the rest has to be binned. Then we
% normalise to the length
%time_start = datetime('now');

theta_i = zeros(n_rays, 1);
phi_i = zeros(n_rays, 1);
theta_o = cell(n_rays, 1);
phi_o = cell(n_rays, 1);

%poolobj = gcp;
%addAttachedFiles(poolobj,'vars.mat')
for i1 = 1:n_rays
    theta=theta_max*rand(1);
    phi=phi_max*rand(1);
    %phi=0;
    theta_i(i1) = theta;
    %length(theta_intervals(theta_intervals<theta));
    %phi_int_loop = phi_intervals{theta_i(i1)};
    phi_i(i1) = phi;
    %length(phi_int_loop(phi_int_loop<phi));
    results = scan_xy_F(tr, usemean, FN, theta, phi, n_x, n_y, n_1, n_2);
    %theta_n = histcounts(results(:,1), theta_intervals);
    theta_o{i1} = results(:,1);
    %phi_n= histcounts(results(:,2), phi_int_loop);
    phi_o{i1} = results(:,2);
    
    
end
% bin outside the loop
matrix_res=zeros(length(angle_vector) + 1, length(angle_vector));
theta_summary=zeros(length(theta_middle)+1, length(theta_middle));
results_per_angle = length(theta_o{1});
% deal with NaNs
for i2 = 1:n_rays
    % find theta_i bin
    theta_i_bin = length(theta_intervals(theta_intervals<theta_i(i2)));
    phi_bins = phi_intervals{theta_i_bin};
    phi_i_bin = length(phi_bins(phi_bins<phi_i(i2)));
    [~ , in_loc] = min(abs(angle_vector(:,1)-theta_i_bin));
    in_loc = in_loc + phi_i_bin - 1;
    for i3 = 1:results_per_angle
        if isnan(theta_o{i2}(i3))
            out_loc = length(angle_vector) + 1;
            theta_o_bin = length(theta_middle) + 1;
        else
            theta_o_bin = length(theta_intervals(theta_intervals<theta_o{i2}(i3)));
            phi_bins = phi_intervals{theta_o_bin};
            phi_o_bin = length(phi_bins(phi_bins<abs(phi_o{i2}(i3))));
            [~ , out_loc] = min(abs(angle_vector(:,1)-theta_o_bin));
            out_loc = out_loc + phi_o_bin -1;
        end
        matrix_res(out_loc, in_loc) = matrix_res(out_loc, in_loc) + 1;
        theta_summary(theta_o_bin, theta_i_bin) = theta_summary(theta_o_bin, theta_i_bin) + 1;
    end
    
end

m_size = size(matrix_res);

for j1=1:m_size(2)
    if sum(matrix_res(:,j1)) ~= 0
        matrix_res(:,j1)=matrix_res(:,j1)/sum(matrix_res(:,j1)); % normalisation to number of rays
    end
end

m_size = size(theta_summary);
for j1=1:m_size(2)
    if sum(theta_summary(:,j1)) ~= 0
        theta_summary(:,j1)=theta_summary(:,j1)/sum(theta_summary(:,j1)); % normalisation to number of rays
    end
end

%datetime('now')-time_start
end

function[angle_res] = scan_xy_F(tr, usemean, FN, theta, phi, n_x, n_y, n_1, n_2)

Lx = max(tr.Points(:,1));
Ly = max(tr.Points(:,2));
h = max(tr.Points(:,3));
r = abs((h+1)/cosd(theta));

r_a_0 = [r*sind(theta)*cosd(phi); r*sind(theta)*sind(phi); r*cosd(theta)];
%x_in = linspace(Lx/100, Lx-(Lx/100), n_x);
%y_in = linspace(Ly/100, Ly-(Ly/100), n_y);
if usemean == true
    x_in = linspace(Lx/100+(Lx/4), (3*Lx/4)-(Lx/100), n_x);
    y_in = linspace(Ly/100+(Ly/4), (3*Ly/4)-(Ly/100), n_y);
    offset = mean(tr.Points(:,3));
else
    x_in = linspace(Lx/100, Lx-(Lx/100), n_x);
    y_in = linspace(Ly/100, Ly-(Ly/100), n_y);
    offset = 0;
end
[X, Y] = meshgrid(x_in, y_in);
theta_res = zeros(numel(X), 1);
phi_res = zeros(numel(X), 1);

for i1 = 1:numel(X)
    x_loop = X(i1);
    y_loop = Y(i1);
    r_a = r_a_0+[x_loop; y_loop; offset];
    r_b = [x_loop; y_loop; offset];
    
    ray_path = r_a;
    d = (r_b-r_a)/norm(r_b-r_a);
    %d_out = d;
    if theta <= 90
        side = 1;
    else
        side = -1;
    end
    n_int = 0;
    intersect = true;
    mirrored = false;
    absorbed = false;
    new_dir = d; % why does this sometimes give an error?
    while intersect == true
        results = single_check(tr, FN, r_a, d, side);
        
        if all(isnan(results.intersn)) && mirrored == false
            
            m_xy = d(2)/d(1);
            found = false;
            if d(1) > 0
                if m_xy*(Lx - r_a(1))+r_a(2) >= 0 && m_xy*(Lx - r_a(1))+r_a(2) <= Ly
                    %disp('right')
                    r_a = [2*Lx-r_a(1); r_a(2); r_a(3)];
                    d = [-d(1); d(2); d(3)]; % mirror
                    found = true;
                end
            elseif m_xy*(0 - r_a(1))+r_a(2) >= 0 && m_xy*(0 - r_a(1))+r_a(2) <= Ly
                %disp('left')
                r_a = [-r_a(1); r_a(2); r_a(3)];
                d = [-d(1); d(2); d(3)]; % mirror
                found = true;
            end
            
            if d(2) > 0 && found == false
                if (Ly - r_a(2))/m_xy + r_a(1) >= 0 && (Ly - r_a(2))/m_xy + r_a(1) <= Lx
                    %disp('top')
                    r_a = [r_a(1); 2*Ly-r_a(2); r_a(3)];
                    d = [d(1); -d(2); d(3)]; % mirror
                else
                    %disp('bottom')
                    r_a = [r_a(1); -r_a(2); r_a(3)];
                    d = [d(1); -d(2); d(3)]; % mirror
                end
                
            end
            mirrored = true;
            
        elseif all(isnan(results.intersn)) && mirrored == true
            intersect = false; % stop
            
        else
            n_int = n_int + 1;
            absorbed = false;
            intersn = results.intersn;
            theta_i = results.theta;
            % look-up table time...
            N = results.N;
            if side == 1
                ni = n_1;
                nj = n_2;
            else
                ni = n_2;
                nj = n_1;
            end
            theta_c = real(asind(nj/ni));
            if theta_i > theta_c % total internal reflection
                new_dir = d - 2*dot(d, N)*N; % reflection
            else
                rnd = rand();
                theta_j=asind(sind(theta_i)*ni/nj);
                %ts=2*ni*cosd(theta_i)/(ni*cosd(theta_i)+nj*cosd(theta_j));
                %tp=2*ni*cosd(theta_i)/(nj*cosd(theta_i)+ni*cosd(theta_j));
                rs=(ni*cosd(theta_i)-nj*cosd(theta_j))/(ni*cosd(theta_i)+nj*cosd(theta_j));
                rp=(ni*cosd(theta_j)-nj*cosd(theta_i))/(nj*cosd(theta_i)+ni*cosd(theta_j));
                R = 0.5*(abs(rs)^2 + abs(rp)^2);
                if rnd <= R
                    new_dir = d - 2*dot(d, N)*N; % reflection
                else
                    %if TRseq(i2) == 'T' & TIR == false
                    %disp('tranmission')
                    tr_par = ((real(n_1)/real(n_2))^side)*(d-dot(d, N)*N);
                    tr_perp = -sqrt(1-norm(tr_par)^2)*N;
                    new_dir = tr_par+tr_perp;
                    side = -side;
                end
            end
            
            ray_path = [ray_path intersn];
            %plot3([intersn(1) intersn(1)+N(1)/2], [intersn(2) intersn(2)+N(2)/2], ...
            %[intersn(3) intersn(3)+N(3)/2], 'g', 'LineWidth', 2)
            if absorbed == false
                d = new_dir;
                d_out = d;
                r_a = intersn+new_dir/1000;
                mirrored = false;
            else
                intersect = false; % stop, ray has been absorbed
                % don't want to calculate the same intersection again!
            end
            
        end
    end
    if absorbed == false
        theta_res(i1) = acosd(d_out(3)/(norm(d_out)^2));
        phi_res(i1) = atand(d_out(2)/d_out(1));
        %ray_path = [ray_path ray_path(:, end) + new_dir];
    else
        theta_res(i1) = NaN;
        phi_res(i1) = NaN;
    end
    %angle_res(i1, 3) = n_int;
    %plot3(ray_path(1,:), ray_path(2,:), ray_path(3,:), 'r','LineWidth', 2)
end
angle_res = [theta_res phi_res];
end
