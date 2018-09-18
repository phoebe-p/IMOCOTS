% absorbed in front surface:
function [v_0, angle_vector] = calc_v_0(matrix_name, options)
% if phi_in = 'uniform', distribute over all phis for specified theta

front_m = load(matrix_name);
fmnames = fieldnames(front_m);
front_m = front_m.(fmnames{1});

theta_in = options.theta_in;
phi_in = options.phi_in;
%options.c_azimuth = 1;
%theta_in = 9;
%phi_in = 'uniform';

if theta_in == 0
    theta_in = 1e-9;
end

if phi_in == 0
    phi_in = 1e-9;
end

n_angle_bins = options.n_angle_bins;
c_azimuth = options.c_azimuth;

sin_a_b = linspace(0, 1, n_angle_bins+1); % number of bins is between 0 and 90
% will have the same number of
% bins between 90 and 180
% phi_max_geometry = 45;
theta_intervals=[asind(sin_a_b)  180-fliplr(asind(sin_a_b(1:end-1)))];
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
%angle_vector = angle_vector(angle_vector(:, 2) <= options.theta_max, :);
theta_i_bin = length(theta_intervals(theta_intervals<theta_in));
phi_bins = phi_intervals{theta_i_bin};
v_in = zeros(length(angle_vector(:,1)), 1);

if strcmp(phi_in, 'uniform')
    I_per_bin = 1/(length(phi_bins)-1);
    [~ , in_loc] = min(abs(angle_vector(:,1)-theta_i_bin));
    in_loc = in_loc:1:(in_loc+length(phi_bins)-2);
    v_in(in_loc, 1) = I_per_bin;
    
else
    phi_i_bin = length(phi_bins(phi_bins<phi_in));
    [~ , in_loc] = min(abs(angle_vector(:,1)-theta_i_bin));
    in_loc = in_loc + phi_i_bin - 1;
    v_in(in_loc, 1) = 1;
end

v_in = v_in(1:length(v_in(:,1))/2, 1);

angle_vector = angle_vector(1:length(angle_vector(:,1))/2, 2:3);

length(v_in);
n_cols = length(front_m(1,:));

front_m = [front_m; zeros(1, n_cols)]; 

v_0 = front_m(:, 3:end)*v_in;
%save('v_0', 'v_0')

end
% last entry in P_inside is absorbed