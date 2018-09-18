function [A_tot, A_front, A_back, R, T] = calculate_absorption(superstrate, stack, substrate, options)

bulk_layer_material =  load(stack(2).layer(1).mat{1});								
bulk_layer_thickness = stack(2).layer(1).thick/1000; 	

if options.parallel
generate_matrices_par(superstrate, stack, substrate, options)
else
    generate_matrices(superstrate, stack, substrate, options)
end


disp('OPTOS absorptance calculation running...');
iter_or_geom = 'iter';	
power_threshold = options.power_threshold;
% addpaths, load material data and matrix names																			%beginning and end of matrix name (see "matrices" folder)
%only used for result file name {front structure, back structure }

%initialization for absorptance variable
wavelength_list = options.wavelengths;
Abs_save = zeros(length(wavelength_list),1);
Abs_front = zeros(length(wavelength_list),1);
Abs_back = zeros(length(wavelength_list),1);
R = zeros(length(wavelength_list),1);
T = zeros(length(wavelength_list),1);

%wavelength loop
for wavelength = wavelength_list
    n_bulk = interp1(bulk_layer_material(:,1),bulk_layer_material(:,2),wavelength/1000); %bulk medium silicon n
    k_bulk = interp1(bulk_layer_material(:,1),bulk_layer_material(:,3),wavelength/1000); %bulk medium silicon k
    alpha = 4*pi*k_bulk/(wavelength/1000);	%calculation of absorption coefficient
    
    [v0,B,C,points] = load_matrices(wavelength, options, stack(1).name, stack(3).name);
    [D_si] = propagation_matrix_D(bulk_layer_thickness,alpha,points);
    
    v_rear = D_si*v0; 																							
    % remaining power fraction vector at the rear side after one pass through bulk
    if sum(v_rear(1:end-1)) < power_threshold
        %disp('a')
        Abs_save(wavelength_list == wavelength) = sum(v0(1:(end-1)));
        Abs_front(wavelength_list == wavelength) = v0(end-1);% calculate absorptance via v0
        Abs_back(wavelength_list == wavelength) = v0(end);
    else
        %disp('b')
        [Abs_save(wavelength_list == wavelength), ...
            Abs_front(wavelength_list == wavelength), ...
            Abs_back(wavelength_list == wavelength), ...
            R(wavelength_list == wavelength), ...
            T(wavelength_list == wavelength)] = ...
            absorb(iter_or_geom,B,C,D_si,v0,points,power_threshold);
    end
end

A_tot = [wavelength_list', Abs_save];
A_front = [wavelength_list', Abs_front];
A_back = [wavelength_list', Abs_back];
R = [wavelength_list', R];
T = [wavelength_list', T];

end