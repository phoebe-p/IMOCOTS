function [Rmat, Tmat, Amat, angle_vector, theta_summary] = RCWA(SC, grating, options, w_index)
% wavelength is already stored in SC.
% need to build the structure first.

% doesn't really make sense to do this on every wavelength loop.
n_angle_bins = options.n_angle_bins;
theta_in_max = options.theta_max;
phi_in_max = options.phi_max;
sin_a_b = linspace(0, 1, n_angle_bins+1); % number of bins is between 0 and 90
c_azimuth = options.c_azimuth;
% will have the same number of
% bins between 90 and 180
%phi_max_geometry = options.phi_max;
%phi_dep = options.phi_dep;
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

%----------%
d_mat = [grating.d21, grating.d22; grating.d31, grating.d32];
% f_mat*d_mat = I (identity matrix); f_mat = inv(d_mat)w
f_mat =  inv(d_mat);
% f_mat = [f21, f31; f22, f32]
f21 = f_mat(1, 1);
f22 = f_mat(2, 1);
f31 = f_mat(1, 2);
f32 = f_mat(2, 2);

angle_vector_in = angle_vector(angle_vector(:, 2) <= theta_in_max, :);
angle_vector_in = angle_vector_in(angle_vector_in(:,3) >= 0 & angle_vector_in(:,3) <= phi_in_max, :);
%if phi_dep == false
%    [~, which_unique] = unique(angle_vector_in(:,2));
%    angle_vector_in = angle_vector_in(which_unique, :);
%end

matrix_res=zeros(length(angle_vector) + 1, length(angle_vector));
theta_summary=zeros(length(theta_middle) + 1, length(theta_middle));
%save('angle_vector_in', 'angle_vector_in')
%length(angle_vector_in)
for j1 = 1:length(angle_vector_in)
   % j1
    theta_in = angle_vector_in(j1, 2); % from 0 to 90 -> bin for theta_in is j1
    phi_in = angle_vector_in(j1, 3);
    [grating, SC, R, T, inc_field] = RCWA_calc(grating, SC, theta_in/(180/pi), phi_in/(180/pi), w_index);
    %matrix_A_s(j1) = SC.RAT.As(w_index);
    %matrix_A_p(j1) = SC.RAT.Ap(w_index);
    % GdCalc.pdf eqn 3.27
    
    % eqn 4.7
    m1 = extractfield(R, 'm1'); % same for R and T
    m2 = extractfield(R, 'm2');
    f2_r = inc_field.f2 + m1*f21 + m2*f22; % = f2_t
    f3_r = inc_field.f3 + m1*f31 + m2*f32; % = f3_t
    f1_r = real(sqrt((cell2mat(grating.pmt(grating.pmt_sup_index))/...
        (SC.wavelengths(w_index))^2) - f2_r.^2 - f3_r.^2));
    % imaginary f1 means no propagating modes
    f1_t = real(-sqrt((cell2mat(grating.pmt(grating.pmt_sub_index))/...
        (SC.wavelengths(w_index))^2) - f2_r.^2 - f3_r.^2));
    
    % calculate angles
    theta_r = num2cell((180/pi)*acos(f1_r./sqrt(f1_r.^2 + f2_r.^2 + f3_r.^2)));
    % measured from normal pointing into superstrate
    theta_t = num2cell((180/pi)*(pi - acos(-f1_t./sqrt(f1_t.^2 + f2_r.^2 + f3_r.^2))));
    % t angles converted to match format of ray tracer
    phi = num2cell((180/pi)*atan(f3_r./f2_r));
    
    % I think there is some ambiguity here about the quadrant phi is in.
    % Should not matter if you have symmetry about x and y axes (i.e. 30
    % degrees from the x axis the same as 30 degrees from the y axis).
    
    % need to "fold back" into symmetry! get angles in range -90 to 90
    % this way. For now, just make bins from -90 to 90.
    
    theta_i_bin = length(theta_intervals(theta_intervals<theta_in));
    phi_bins_R = phi_intervals{theta_i_bin};
    phi_i_bin = length(phi_bins_R(phi_bins_R<phi_in));
    [~ , in_loc] = min(abs(angle_vector(:,1)-theta_i_bin));
    in_loc = in_loc + phi_i_bin - 1;
    
    for i1 = 1:length(theta_r)
        if isnan(theta_r{i1})
            theta_r{i1} = 0;
        end
        if isnan(theta_t{i1})
            theta_t{i1} = 0;
        end
        if phi{i1} < 0
            phi{i1} = phi{i1} + 90;
        end
        
        theta_o_bin_R = length(theta_intervals(theta_intervals<theta_r{i1}));
        phi_bins_R = phi_intervals{theta_o_bin_R};
        phi_o_bin_R = length(phi_bins_R(phi_bins_R<abs(phi{i1})));
        [~ , out_loc_R] = min(abs(angle_vector(:,1)-theta_o_bin_R));
        out_loc_R = out_loc_R + phi_o_bin_R - 1;
        R_eff = 0.5*(R(i1).eff1 + R(i1).eff2); % UNPOLARIZED
        matrix_res(out_loc_R, in_loc) = matrix_res(out_loc_R, in_loc) + R_eff;
        theta_summary(theta_o_bin_R, theta_i_bin) = theta_summary(theta_o_bin_R, theta_i_bin) + R_eff;
        
        theta_o_bin_T = length(theta_intervals(theta_intervals<theta_t{i1}));
        phi_bins_T = phi_intervals{theta_o_bin_T};
        phi_o_bin_T = length(phi_bins_T(phi_bins_T<abs(phi{i1})));
        [~ , out_loc_T] = min(abs(angle_vector(:,1)-theta_o_bin_T));
        out_loc_T = out_loc_T + phi_o_bin_T - 1;
        T_eff = 0.5*(T(i1).eff1 + T(i1).eff2); % UNPOLARIZED
        matrix_res(out_loc_T, in_loc) = matrix_res(out_loc_T, in_loc) + T_eff;
        theta_summary(theta_o_bin_T, theta_i_bin) = theta_summary(theta_o_bin_T, theta_i_bin) + T_eff;
       
    end
    %0.5*sum(extractfield(R, 'eff1') + extractfield(R, 'eff2'))+ ...
     %        0.5*sum(extractfield(T, 'eff1') + extractfield(T, 'eff2') )
    
end
Rmat = matrix_res(1:(length(angle_vector)/2), 1:(length(angle_vector)/2));
Tmat = flipud(matrix_res((length(angle_vector)/2 + 1):(end-1), 1:(length(angle_vector)/2)));

Amat = zeros(1, length(angle_vector)/2);
% Amat = 1 - sum(Rmat, 1) - sum(Tmat, 1); % this doesn't work if not all
% the 'in' angles are filled, will get erroneous unity absorption!

for i1 = 1:(length(angle_vector)/2)
    if sum(Rmat(:, i1) + Tmat(:, i1)) > 0
        Amat(1, i1) = 1 - sum(Rmat(:, i1) + Tmat(:, i1));
    end
end

m_size = size(theta_summary);
for j1=1:m_size(2)
    if sum(theta_summary(:,j1)) ~= 0
        theta_summary(:,j1)=theta_summary(:,j1)/sum(theta_summary(:,j1)); % normalisation to number of rays
    end
end
% theta_summary DOES NOT make sense! is wrong!

% not actually binning absorption in matrix_res!
%Amat = matrix_res(end, 1:length(angle_vector)/2);
end


function [grating, SC, R, T, inc_field] = RCWA_calc(grating, SC, theta_in, phi_in, index)

SC.theta=theta_in; % incidence polar angle from x1 axis
SC.phi=phi_in;
%plot_grating(grating, SC.period)
[grating, SC, R, T, inc_field]=optical(grating, SC, index);


end

function [grating, SC, R, T, inc_field]=optical(grating, SC, index)
j1 = index;
layer_indices = needed_layers(SC);
grating.pmt={};
grating.pmt{1}=SC.layer(1).material.pmt(j1);
    for j2=2:length(SC.layer)
        for j3=1:length(SC.layer(j2).material)
            grating.pmt{end+1}=SC.layer(j2).material(j3).pmt(j1);
        end        
    end
    %calculate scattered fields - internal and external
    [grating, SC, R, T, inc_field] = scatfield(SC, grating, j1);
    [grating, SC] = Poynting_abs(SC, grating, j1, layer_indices);
end

function[grating, SC] = Poynting_abs(SC, grating, index, layer_indices)
        calc_layer = layer_indices.calc_layer;
        t = layer_indices.t;
        b = layer_indices.b;
        
        P_needed = sort(unique([b'; t']));
        P_data = table(calc_layer', t', b', 'VariableNames',{'i', 't', 'b'});
        for j = 1:length(P_needed)
            a = struct2table(grating.stratum{1, P_needed(j)}.full_field);
            ffE2s = table2array(a(:, 7));
            ffE3s = table2array(a(:, 8));
            ffH2s = table2array(a(:, 10));
            ffH3s = table2array(a(:, 11));
            
            ffE2p = table2array(a(:, 13));
            ffE3p = table2array(a(:, 14));
            ffH2p = table2array(a(:, 16));
            ffH3p = table2array(a(:, 17));

            P_needed(j, 2) = sum(real(ffE2s.*conj(ffH3s) - ffE3s.*conj(ffH2s)));
            P_needed(j, 3) = sum(real(ffE2p.*conj(ffH3p) - ffE3p.*conj(ffH2p)));
            % p polarization!
            
            %zeroth order:
            % zeroth order; 
            ind = find(table2array(a(:,1)) == 0 & table2array(a(:,2)) == 0);
            P_needed(j, 4) = real(ffE2s(ind).*conj(ffH3s(ind)) - ffE3s(ind).*conj(ffH2s(ind)));
            P_needed(j, 5) = real(ffE2p(ind).*conj(ffH3p(ind)) - ffE3p(ind).*conj(ffH2p(ind)));
            %%%%----%
            
        end
        findP_s = @(x) - P_needed(P_needed(:,1) == x, 2);
        findP_p = @(x) - P_needed(P_needed(:,1) == x, 3);
        findP_s_zeroth = @(x) - P_needed(P_needed(:,1) == x, 4);
        findP_p_zeroth = @(x) - P_needed(P_needed(:,1) == x, 5);
        %%

        P_data.P_abs_s = (arrayfun(findP_s, P_data.t) - arrayfun(findP_s, P_data.b))/cos(SC.theta);
        P_data.P_abs_p = (arrayfun(findP_p, P_data.t) - arrayfun(findP_p, P_data.b))/cos(SC.theta);
        P_data.P_abs_s_zeroth = (arrayfun(findP_s_zeroth, P_data.t) - arrayfun(findP_s_zeroth, P_data.b))/cos(SC.theta);
        P_data.P_abs_p_zeroth = (arrayfun(findP_p_zeroth, P_data.t) - arrayfun(findP_p_zeroth, P_data.b))/cos(SC.theta);
        for j2 = 1:length(P_data.i)
        SC.layer(P_data.i(j2)).P_abs_s(index) = P_data.P_abs_s(j2);
        SC.layer(P_data.i(j2)).P_abs_p(index) = P_data.P_abs_p(j2);
        SC.layer(P_data.i(j2)).P_abs(index) = (P_data.P_abs_s(j2) + P_data.P_abs_p(j2))/2;
        
        SC.layer(P_data.i(j2)).P_abs_s_zeroth(index) = P_data.P_abs_s_zeroth(j2);
        SC.layer(P_data.i(j2)).P_abs_p_zeroth(index) = P_data.P_abs_p_zeroth(j2);
        SC.layer(P_data.i(j2)).P_abs_zeroth(index) = (P_data.P_abs_s_zeroth(j2) + P_data.P_abs_p_zeroth(j2))/2;
        end
  
end

%----------------%
