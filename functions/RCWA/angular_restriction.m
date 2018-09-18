function SC = angular_restriction(grating, SC, R_w, T_w, inc_field)
%[grating, SC, R_w, T_w, inc_field] = square_grid_2(173, 179, 5);

%----------%
d_mat = [grating.d21, grating.d22; grating.d31, grating.d32];
% f_mat*d_mat = I (identity matrix); f_mat = inv(d_mat)
f_mat =  inv(d_mat);
% f_mat = [f21, f31; f22, f32]
f21 = f_mat(1, 1);
f22 = f_mat(2, 1);
f31 = f_mat(1, 2);
f32 = f_mat(2, 2);
max_angle = asin(0.45); % radians

for w_index = 1:length(SC.wavelengths)
    
    R = squeeze(cell2mat(struct2cell(R_w{w_index})))';
    T = squeeze(cell2mat(struct2cell(T_w{w_index})))';
    % eqn 4.7
    m1 = R(:, 1); % same for R and T
    m2 = R(:, 2);
    f2_r = inc_field.f2 + m1*f21 + m2*f22; % = f2_t
    f3_r = inc_field.f3 + m1*f31 + m2*f32; % = f3_t
    f1_r = real(sqrt((cell2mat(grating.pmt(grating.pmt_sup_index))/...
        (SC.wavelengths(w_index))^2) - f2_r.^2 - f3_r.^2));
    % imaginary f1 means no propagating modes
    f1_t = real(-sqrt((cell2mat(grating.pmt(grating.pmt_sub_index))/...
        (SC.wavelengths(w_index))^2) - f2_r.^2 - f3_r.^2));
    
    % calculate angles
    theta_r = acos(f1_r./sqrt(f1_r.^2 + f2_r.^2 + f3_r.^2));
    % measured from normal pointing into superstrate
    theta_t = acos(-f1_t./sqrt(f1_t.^2 + f2_r.^2 + f3_r.^2));
    % measured from normal pointing into substrate
    
    theta_r(isnan(theta_r)) = 0;
    theta_t(isnan(theta_t)) = 0;
    
    R_zeroth = R(R(:,1) == 0 & R(:,2) == 0, :);
    T_zeroth = T(T(:,1) == 0 & T(:,2) == 0, :);
    R_res = R(theta_r < max_angle, :);
    T_res = T(theta_r < max_angle, :);

    %store the total reflection, absorption and transmission for s an p
    %incidence. The tranmsission is into the bottom cell. The absorption is
    %everywhere else (middle and top cells, window, ARC and NPs)
    SC.RAT.Rs_res(w_index)=sum(R_res(:, 3));
    SC.RAT.Rp_res(w_index)=sum(R_res(:, 4));
    SC.RAT.Rtotal_res(w_index)=sum((R_res(:,3)+R_res(:,4))/2);
    SC.RAT.Ts_res(w_index)=sum(T_res(:, 3));
    SC.RAT.Tp_res(w_index)=sum(T_res(:, 4));
    SC.RAT.Ttotal_res(w_index)=sum((T_res(:, 3)+T_res(:, 4))/2);
    
    SC.RAT.Rs_zeroth(w_index)=sum(R_zeroth(:, 3));
    SC.RAT.Rp_zeroth(w_index)=sum(R_zeroth(:, 4));
    SC.RAT.Rtotal_zeroth(w_index)=sum((R_zeroth(:,3)+R_zeroth(:,4))/2);
    SC.RAT.Ts_zeroth(w_index)=sum(T_zeroth(:, 3));
    SC.RAT.Tp_zeroth(w_index)=sum(T_zeroth(:, 4));
    SC.RAT.Ttotal_zeroth(w_index)=sum((T_zeroth(:, 3)+T_zeroth(:, 4))/2);
    
end

%figure
%plot(SC.wavelengths, SC.RAT.Rtotal, '-', SC.wavelengths, SC.RAT.Rtotal_res, '--', ...
%    SC.wavelengths, SC.RAT.Rtotal_zeroth, '.')
