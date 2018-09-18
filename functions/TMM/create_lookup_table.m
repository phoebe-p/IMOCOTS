function [] = create_lookup_table(superstrate, layers, substrate, options, name)
% create lookup table

% consider the effect of non-real refractive index (absorption) in
% incidence and transmission medium?
absorbing = options.incident_absorption;

addpath('nk', 'functions')
theta_in = (0:options.lookuptable_step:options.theta_max)*pi/180;
%wavelengths = min(options.wavelengths):10:max(options.wavelengths);
wavelengths = options.wavelengths;

% if the incidence medium is very highly absorbing (like a metal), get very
% weird results. Set a threshold: if alpha (um-1) > 100, just assume that
% any ray which is already in that material will stay there (i.e. will
% eventually get absorbed there; so R = 1.

incidence_medium.material_filename=superstrate;
incidence_medium.refractive_index=[];
%
stack = [];

if isfield(layers, 'eff_mat') == 0
    for i1 = 1:(length(layers))
        layers(i1).eff_mat = 0;
    end
end
for i1 = 1:(length(layers))
    if isempty(layers(i1).eff_mat)
        layers(i1).eff_mat = 0;
    end
end
for i1 = 1:(length(layers))
    if layers(i1).eff_mat == 0
        stack(end+1).material_filename=layers(i1).mat{1};
        stack(end).refractive_index=[];
        stack(end).thickness=layers(i1).thickness;
    else
        stack = make_effective_medium(stack, layers(i1).mat{1}, layers(i1).mat{2}, layers(i1).n_layers, layers(i1).thickness, wavelengths);
    end
end

behind_medium.material_filename=substrate;
behind_medium.refractive_index=[];

[incidence_medium]=get_refractive_index_TMM(incidence_medium, wavelengths);
[behind_medium]=get_refractive_index_TMM(behind_medium, wavelengths);

if ~absorbing
    incidence_medium.refractive_index = real(incidence_medium.refractive_index);
    behind_medium.refractive_index = real(behind_medium.refractive_index);
end

n1 = incidence_medium.refractive_index;
n2 = behind_medium.refractive_index;

% loop for light incident from front
[R_plus, T_plus, A_plus, R_minus, T_minus, A_minus, alpha_plus, alpha_minus] = ...
    deal(zeros(length(theta_in), length(wavelengths)));

for i1 = 1:length(theta_in)
    [R_loop_s, A_loop_s, T_loop_s, R_loop_p, A_loop_p, T_loop_p, alpha_out] = ...
        TMM_lookuptable(wavelengths, theta_in(i1), incidence_medium, stack, behind_medium);
    R_plus(i1, :) = (R_loop_s+R_loop_p)/2;
    T_plus(i1, :) = (T_loop_s+T_loop_p)/2;
    A_plus(i1, :) = (A_loop_s+A_loop_p)/2;
    alpha_plus(i1, :) = alpha_out;
end

b = incidence_medium;
incidence_medium = behind_medium;
behind_medium = b;
stack = flip(stack);

for i1 = 1:length(theta_in)
    [R_loop_s, A_loop_s, T_loop_s, R_loop_p, A_loop_p, T_loop_p, alpha_out] = ...
        TMM_lookuptable(wavelengths, theta_in(i1), incidence_medium, stack, behind_medium);
    R_minus(i1, :) = (R_loop_s+R_loop_p)/2;
    T_minus(i1, :) = (T_loop_s+T_loop_p)/2;
    A_minus(i1, :) = (A_loop_s+A_loop_p)/2;
    alpha_minus(i1, :) = alpha_out;
end

% fi label in filename: Front Incidence
save(strcat('results/', options.name, '/', name, '_fi'), 'R_plus', 'T_plus', ...
    'A_plus', 'R_minus', 'T_minus', 'A_minus', 'theta_in', 'n1', 'n2', ...
    'wavelengths', 'alpha_plus', 'alpha_minus')

% ri label in filename: Rear Incidence
a1 = R_minus;
a2 = T_minus;
a3 = A_minus;

R_minus = R_plus;
T_minus = T_plus;
A_minus = A_plus;

R_plus = a1;
T_plus = a2;
A_plus = a3;
save(strcat('results/', options.name, '/', name, '_ri'), 'R_plus', 'T_plus', ...
    'A_plus', 'R_minus', 'T_minus', 'A_minus', 'theta_in', 'n1', 'n2', ...
    'wavelengths', 'alpha_plus', 'alpha_minus')
end


function [Rs, As, Ts, Rp, Ap, Tp, alpha_o, stack]=TMM_lookuptable(wavelengths, theta, incidence_medium, stack, behind_medium)

%do not touch from this point on ------------------------------------------
th_r=theta;
%[incidence_medium]=get_refractive_index_TMM(incidence_medium, wavelengths);
[stack]=get_refractive_index_TMM(stack, wavelengths);
%save('stack', 'stack')
%[behind_medium]=get_refractive_index_TMM(behind_medium, wavelengths);
Ts = zeros(length(wavelengths), 1);
Tp = zeros(length(wavelengths), 1);
Rs = zeros(length(wavelengths), 1);
Rp = zeros(length(wavelengths), 1);
alpha_o = zeros(length(wavelengths), 1);

A_layer_s = zeros(length(wavelengths), length(stack));
A_layer_p = zeros(length(wavelengths), length(stack));

P_enter_s = zeros(length(wavelengths), 1);
P_enter_p = zeros(length(wavelengths), 1);

for j1=1:length(wavelengths)
    threshold = 1e-5;
    wavelength=wavelengths(j1);   %vacuum wavelength (nm)
    if 4*pi*imag(incidence_medium.refractive_index(j1))/wavelength > 0.01
        % highly absorbing incidence medium
        %disp(['Lambda =  ' num2str(wavelength) ' - highly absorbing incidence medium, set R = 1'])
        alpha_o(j1) = pi/2;
        Ts(j1)=0;
        Tp(j1)=0;
        P_enter_s = 0; % so that R = 1
        P_enter_p = 0;
        A_layer_s(j1, :)=0;
        A_layer_p(j1, :)=0;
    else
        theta_incidence_rad = complex_incidence(th_r, incidence_medium.refractive_index(j1));
        
        if imag(incidence_medium.refractive_index(j1)*cos(theta_incidence_rad)) < -threshold
            % i.e. real part of cos(theta_j) is basically zero
            theta_incidence_rad = pi-theta_incidence_rad;
        end
        
        [R_s, R_p, T_s, T_p, layer_matrices, t_s, t_p, r_s, r_p, th_out]=stack_RAT(incidence_medium, stack, behind_medium, theta_incidence_rad, wavelength, j1);
        alpha_o(j1) = th_out;
        Ts(j1)=T_s;
        Tp(j1)=T_p;
        Rs(j1)=R_s;
        Rp(j1)=R_p;
        %As(j1)=1-Ts(j1)-Rs(j1);
        %Ap(j1)=1-Tp(j1)-Rp(j1);
        
        ni=incidence_medium.refractive_index(j1);
        
        [s_abs, p_abs] = layer_absorptions(layer_matrices, t_s, t_p, stack, wavelength, j1, ni, theta_incidence_rad);
        A_layer_s(j1, :)=s_abs;
        A_layer_p(j1, :)=p_abs;
        
        P_enter_s(j1) = 1-Rs(j1) + 2*imag(r_s)*imag(ni*cos(theta_incidence_rad))/real(ni*cos(theta_incidence_rad));
        P_enter_p(j1) = 1-Rp(j1) - 2*imag(r_p)*imag(ni*cos(conj(theta_incidence_rad)))/real(ni*cos(conj(theta_incidence_rad)));
    end
    % if any([R_s, R_p, T_s, T_p] < -1e-8) || any([R_s, R_p, T_s, T_p] > 1 + 1e-8)
    %end
end

As = sum(A_layer_s, 2);
Ap = sum(A_layer_p, 2);
Rs = 1-P_enter_s;
Rp = 1-P_enter_p;
%T=(Ts+Tp)/2;
%R =(Rs+Rp)/2;
%A=(As+Ap)/2;
%P_enter = (P_enter_s + P_enter_p)/2; % instead of R, use 1-P_enter so that everything always sums to 1.
%A_layer = (A_layer_s + A_layer_p)/2;
% sumAs = sum(A_layer_s, 2);
% sumAp = sum(A_layer_p, 2);
% effRs = 1 - P_enter_s;
% effRp = 1 - P_enter_p;

end

function [R_s, R_p, T_s, T_p, layer_matrices, t_s, t_p, r_s, r_p, theta_out]=stack_RAT(previous_medium, stack, next_medium, theta_incidence, wavelength, j1)
threshold = 1e-5;

n0=previous_medium.refractive_index(j1);
theta_0=theta_incidence;
%stack_matrix_s=eye(2);
%stack_matrix_p=eye(2);
% n = 0; incidence material
theta_n = theta_incidence;
n_n = previous_medium.refractive_index(j1);
n_nplus1 = stack(1).refractive_index(j1);
%ni=previous_medium.refractive_index(j1);
%    theta_i=theta_incidence;

theta_nplus1 = asin(sin(theta_n)*n_n/n_nplus1);

% if abs(imag(theta_nplus1)) < threshold % this is to avoid weird errors with
%     % very small imaginary components which sometimes cause theta_nplus1 to
%     % oscillate and give spurious results
%     theta_nplus1 = real(theta_nplus1);
% end
if imag(n_nplus1*cos(theta_nplus1)) < -threshold
    % i.e. real part of cos(theta_j) is basically zero
    theta_nplus1 = pi-theta_nplus1;
end
ts=2*n_n*cos(theta_n)/(n_n*cos(theta_n)+n_nplus1*cos(theta_nplus1));
tp=2*n_n*cos(theta_n)/(n_nplus1*cos(theta_n)+n_n*cos(theta_nplus1));
rs=(n_n*cos(theta_n)-n_nplus1*cos(theta_nplus1))/(n_n*cos(theta_n)+n_nplus1*cos(theta_nplus1));
rp=(n_nplus1*cos(theta_n)-n_n*cos(theta_nplus1))/(n_nplus1*cos(theta_n)+n_n*cos(theta_nplus1));
stack_matrix_s = [1 rs; rs 1]*1/ts;
stack_matrix_p = [1 rp; rp 1]*1/tp;

theta_n = theta_nplus1;
n_n = n_nplus1;

for j2=1:length(stack)
    if j2==length(stack)
        n_nplus1=next_medium.refractive_index(j1);
    else
        n_nplus1=stack(j2+1).refractive_index(j1);
    end
    
    theta_nplus1=asin(sin(theta_n)*n_n/n_nplus1);
    %     if abs(imag(theta_nplus1)) < threshold
    %         theta_nplus1 = real(theta_nplus1);
    %     end
    if imag(n_nplus1*cos(theta_nplus1)) < -threshold
        % i.e. real part of cos(theta_j) is basically zero
        theta_nplus1 = pi-theta_nplus1;
    end
    % calculate coefficients for n -> n+1
    ts=2*n_n*cos(theta_n)/(n_n*cos(theta_n)+n_nplus1*cos(theta_nplus1));
    tp=2*n_n*cos(theta_n)/(n_nplus1*cos(theta_n)+n_n*cos(theta_nplus1));
    rs=(n_n*cos(theta_n)-n_nplus1*cos(theta_nplus1))/(n_n*cos(theta_n)+n_nplus1*cos(theta_nplus1));
    rp=(n_nplus1*cos(theta_n)-n_n*cos(theta_nplus1))/(n_nplus1*cos(theta_n)+n_n*cos(theta_nplus1));
    delta=2*pi*stack(j2).thickness*n_n*cos(theta_n)/wavelength;
    %if imag(delta) > 35:
    %       delta = real(delta) + 35*1i
    %end
    layer_matrices(j2).Ms_interface=[1,rs;rs,1]/ts;
    layer_matrices(j2).Mp_interface=[1,rp;rp,1]/tp;
    layer_matrices(j2).M_interior=[exp(-1i*delta),0;0,exp(1i*delta)];
    layer_matrices(j2).theta=theta_n; % ANGLE IN
    
    stack_matrix_s=stack_matrix_s*layer_matrices(j2).M_interior*layer_matrices(j2).Ms_interface;
    stack_matrix_p=stack_matrix_p*layer_matrices(j2).M_interior*layer_matrices(j2).Mp_interface;
    
    theta_n = theta_nplus1;
    n_n = n_nplus1;
end

%stack_intensity_matrix_s=abs(stack_matrix_s).*abs(stack_matrix_s)*real(previous_medium.refractive_index(j1)*cos(theta_incidence))/real(next_medium.refractive_index(j1)*cos(theta_next));
%stack_intensity_matrix_p=abs(stack_matrix_p).*abs(stack_matrix_p)*real(previous_medium.refractive_index(j1)*cos(theta_incidence))/real(next_medium.refractive_index(j1)*cos(theta_next));

t_s = 1/stack_matrix_s(1,1);
t_p = 1/stack_matrix_p(1,1);

r_s = stack_matrix_s(2,1)/stack_matrix_s(1,1);
r_p = stack_matrix_p(2,1)/stack_matrix_p(1,1);

R_s = abs(r_s)^2;
R_p = abs(r_p)^2;

T_s = abs(t_s)^2*real(n_nplus1*cos(theta_nplus1))/real(n0*cos(theta_0));
T_p = abs(t_p)^2*real(n_nplus1*cos(conj(theta_nplus1)))/real(n0*cos(conj(theta_0)));
%theta_out = theta_nplus1;
% we want to know the angle at which the real Poynting vector is flowing
% out; tan(alpha) = k_x/real(k_z). If n_0 and n_nplus1 are real, this is
% just theta_nplus1
theta_out = atan(real(n_nplus1*sin(theta_nplus1))/real(n_nplus1*cos(theta_nplus1)));

end

function [medium]=get_refractive_index(medium, wavelengths)
for j1=1:length(medium)
    if isempty(medium(j1).material_filename)
        medium(j1).refractive_index=medium(j1).refractive_index*ones(1,length(wavelengths));
    else
        nk=load(medium(j1).material_filename);
        nk(:,1)=nk(:,1)*1000; %convert from um to nm
        medium(j1).refractive_index=interp1(nk(:,1),nk(:,2),wavelengths)+1i*interp1(nk(:,1),nk(:,3),wavelengths);
    end
end

end

function [s_abs, p_abs] = layer_absorptions(l_mat, t_s, t_p, stack, lambda, wi, ni, theta_i)
% calculate v and w
s_abs = zeros(1, length(stack));
p_abs = zeros(1, length(stack));

vw_s = [t_s; 0];
vw_p = [t_p; 0];


for j1=length(stack):-1:1
    n = stack(j1).refractive_index(wi);
    theta = l_mat(j1).theta;
    k_z = 2*pi*n*cos(theta)/lambda;
    normf = real(ni*cos(theta_i));
    
    Mn_s = l_mat(j1).M_interior*l_mat(j1).Ms_interface;
    vw_s = Mn_s*vw_s;
    l_mat(j1).vw_s = vw_s;
    
    A1_s = imag(n*cos(theta)*k_z)*abs(vw_s(2))^2/normf;
    A2_s = imag(n*cos(theta)*k_z)*abs(vw_s(1))^2/normf;
    A3_s = imag(n*cos(theta)*k_z)*vw_s(1)*conj(vw_s(2))/normf;
    
    %depths = linspace(0, stack(j1).thickness, 100);
    %a_z = zeros(length(depths), 1);
    %for j2 = 1:length(depths)
    %    z = depths(j2);
    %    a_z(j2) = depth_dep_absorption(z, A1_s, A2_s, A3_s, k_z);
    %end
    %plot(depths, a_z)
    %hold on
    
    A_L = @(z) depth_dep_absorption(z, A1_s, A2_s, A3_s, k_z);  % a function handle
    
    s_abs(j1) = integral(A_L, 0, stack(j1).thickness);
    
    Mn_p = l_mat(j1).M_interior*l_mat(j1).Mp_interface;
    vw_p = Mn_p*vw_p;
    l_mat(j1).vw_p = vw_p;
    
    normf = real(ni*cos(conj(theta_i)));
    
    A1_p = 2*imag(k_z)*real(n*cos(conj(theta)))*abs(vw_p(2))^2/normf;
    A2_p = 2*imag(k_z)*real(n*cos(conj(theta)))*abs(vw_p(1))^2/normf;
    A3_p = -2*real(k_z)*imag(n*cos(conj(theta)))*vw_p(1)*conj(vw_p(2))/normf;
    % There is a - sign missing in Steve Byrnes' pdf for A3_p!!!
    
    %a_z = zeros(length(depths), 1);
    %for j2 = 1:length(depths)
    %    z = depths(j2);
    %    a_z(j2) = depth_dep_absorption(z, A1_p, A2_p, A3_p, k_z);
    %end
    %plot(depths, a_z)
    
    A_L = @(z) depth_dep_absorption(z, A1_p, A2_p, A3_p, k_z);  % a function handle
    
    p_abs(j1) = integral(A_L, 0, stack(j1).thickness);
end

end

function [a_z] = depth_dep_absorption(z, A1, A2, A3, k_z)
a_z = A1*exp(2*z*imag(k_z)) + A2*exp(-2*z*imag(k_z)) + ...
    A3*exp(2*1i*z*real(k_z)) + conj(A3)*exp(-2*1i*z*real(k_z));
end


function [theta] = complex_incidence(alpha, ntilde)
n = real(ntilde);
k = imag(ntilde);
D = k^2-n^2;
T = 1+tan(alpha)^2;
% take positive since root1 > -D
kzr = sqrt((-D + sqrt(D^2+4*T*n^2*k^2))/(2*T));
kzi = n*k/kzr;
% n*cos(theta) = kz
theta = acos((kzr+1i*kzi)/ntilde);
end