function [] = generate_matrices(superstrate, stack, substrate, options)

% should be able to make angle_vector beforehand, binning should be the
% same for all

mkdir(strcat('results/', options.name))
addpath(strcat('results/', options.name))
lambda = options.wavelengths;
%warning('off')

% stack(2).type should == 'bulk'.
% argument order: incidence medium, layers, transmission medium

if ~strcmp(stack(2).type, 'bulk')
    error('stack(2).type should be "bulk"');
else
    bulk = stack(2).layer(1).mat{1};
end

% if the front surface is to be treated with ray tracing and requires TMM
% to calculate R, A, T for the layers, need to consider both sides; to
% calculate the redistribution matrix at the front surface for incident
% light and subsequent internal reflections from the bulk.
if strcmp(stack(1).type, 'ray')
    if ~isempty(stack(1).layer)
        name = stack(1).name;
        if exist(strcat(options.name, '/', name, '_fi.mat')) == 0 && ...
                exist(strcat(options.name, '/', name, '_ri.mat')) == 0
            disp('Creating lookup tables for front surface...')
            create_lookup_table(superstrate, stack(1).layer, bulk, options, name)
        else
            disp('Found existing lookup tables for front surface...')
        end
    end
    if ischar(stack(1).texture)
        [stack(1).texture, stack(1).inv_texture, options] = make_texture(stack(1).texture, options);
    end
end

% if the back surface requires TMM, only care about internal reflection (no
% light incident from the bottom of the cell)
if strcmp(stack(3).type, 'ray')
    if ~isempty(stack(3).layer)
        name = stack(3).name;
        if exist(strcat(options.name, '/', name, '_fi.mat')) == 0 && ...
                exist(strcat(options.name, '/', name, '_ri.mat')) == 0
            disp('Creating lookup tables for back surface...')
            create_lookup_table(bulk, stack(3).layer, substrate, options, name)
        else
            disp('Found existing lookup tables for back surface...')
        end
    end
    if ischar(stack(3).texture)
        [stack(3).texture, stack(3).inv_texture, options] = make_texture(stack(3).texture, options);
    end
end


% if the front surface is to be treated with RCWA, need to consider both
% sides; to calculate the redistribution matrix at the front surface for incident
% light and subsequent internal reflections from the bulk
% pre-generate grating and SC structures for feeding into GD-Calc code.
if strcmp(stack(1).type, 'RCWA')
    disp('Creating GD-Calc structures for front surface...')
    [grating_1_1, SC_1_1] = make_grating(superstrate, stack(1), bulk, options);
    flip_layer = flip(stack(1).layer);
    flip_stack = stack(1);
    flip_stack.layer = flip_layer;
    % Need to make sure this works for layers that are not homogeneous in z
    % direction (e.g. pyramids)
    
    [grating_1_2, SC_1_2] = make_grating(bulk, flip_stack, superstrate, options);
end

% if the back surface requires RCWA, only care about internal reflection (no
% light incident from the bottom of the cell)
if strcmp(stack(3).type, 'RCWA')
    disp('Creating GD-Calc structures for back surface...')
    [grating_3_1, SC_3_1] = make_grating(bulk, stack(3), substrate, options);
end


disp('Begin calculating matrices at each wavelength...')
for i1 = 1:length(lambda)
    disp(['Wavelength ' num2str(i1) ' of ', num2str(length(lambda)), ':'])
    if exist(strcat(options.name, '/A_', stack(1).name, num2str(lambda(i1)), 'nm_sum.mat')) == 0

    % front surface:  need to calculate A for incident light at
    % the required angle, and then B for redistribution at the front
    % surface for light incident from below (the bulk)
    if strcmp(stack(1).type, 'ray')
        
        if isempty(stack(1).layer)
            % for incidence from front of structure:
            % A_R: reflected from front surface
            % A_T: transmitted through front surface
            % A_A: absorbed by front surface
            disp('   Front surface Fresnel ray-tracing')
            [A_R, A_T, A_A, angle_vector, A_summary] = RT_Fresnel(superstrate, bulk, stack(1).texture, options, lambda(i1));
            
            % for incidence from inside of structure:
            % B_R: reflected from front surface
            % B_T: transmitted through front surface
            % B_A: absorbed by front surface
            [B_R, B_T, B_A, ~, B_summary] = RT_Fresnel(bulk, superstrate, stack(1).inv_texture, options, lambda(i1));
            % Fresnel calculations will always give zero absorption
            
        else
            disp('   Front surface lookup table (TMM) ray-tracing')
            [A_R, A_T, A_A, angle_vector, A_summary] = RT_Table(stack(1).texture, options, lambda(i1), strcat(stack(1).name, '_fi'));
            [B_R, B_T, B_A,  ~, B_summary] = RT_Table(stack(1).inv_texture, options, lambda(i1), strcat(stack(1).name, '_ri'));
        end
        
    elseif strcmp(stack(1).type, 'RCWA')
        disp('   Front surface RCWA')
        [A_R, A_T, A_A, angle_vector, A_summary] = RCWA(SC_1_1, grating_1_1, options, i1);
        [B_R, B_T, B_A, ~, B_summary] = RCWA(SC_1_2, grating_1_2, options, i1);
        
    elseif strcmp(stack(1).type, 'TMM')    
       disp('   Front surface TMM (planar surface)')
        [A_R, A_T, A_A, angle_vector, A_summary] = TMM(superstrate, stack(1).layer, bulk, options, lambda(i1));
        [B_R, B_T, B_A, ~, B_summary] = TMM(bulk, flip(stack(1).layer), superstrate, options, lambda(i1));
        
    end
    disp('   Done, saving front matrices')
    %save('full_angle_vector.mat', 'angle_vector')
    angle_vector = angle_vector(1:length(angle_vector(:,1))/2, 2:3);
    A = [angle_vector A_T; ...
        0 0 A_A];
    % last row of A = absorption of light (arbitrarily set angle 0)
    mat = A;
    %save(strcat(options.name, '/A_', num2str(lambda(i1)), 'nm.mat'), 'A');
    parsave(strcat(options.name, '/A_', stack(1).name, num2str(lambda(i1)), 'nm.mat'), mat);
    %save(strcat(options.name, '/A_', num2str(lambda(i1)), 'nm_sum.mat'), 'A_summary');
    mat= A_summary;
    parsave(strcat(options.name, '/A_', stack(1).name, num2str(lambda(i1)), 'nm_sum.mat'), mat);
    
    R_front = A_R;
    B = [angle_vector B_R; ...
        0 0 B_A];
    mat = B;
    %save(strcat(options.name, '/B_', num2str(lambda(i1)), 'nm.mat'), 'B');
    parsave(strcat(options.name, '/B_', stack(1).name, num2str(lambda(i1)), 'nm.mat'), mat);
    %save(strcat(options.name, '/B_', num2str(lambda(i1)), 'nm_sum.mat'), 'B_summary');
    % back surface:  need to calculate C for light incident from the bulk
    mat = B_summary;
    parsave(strcat(options.name, '/B_', stack(1).name, num2str(lambda(i1)), 'nm_sum.mat'), mat);
    else disp('Front matrix already calculated')
    end
    
    
    if exist(strcat(options.name, '/C_', stack(3).name, num2str(lambda(i1)), 'nm_sum.mat')) == 0
    if strcmp(stack(3).type, 'ray')
        
        if isempty(stack(3).layer)
            disp('   Back surface Fresnel ray-tracing')
            [C_R, C_T, C_A, angle_vector, C_summary] = RT_Fresnel(bulk, substrate, stack(3).texture, options, lambda(i1));
            % should be able to work out power lost from sum of R and T
        else
            disp('   Back surface lookup table (TMM) ray-tracing')
            [C_R, C_T, C_A, angle_vector, C_summary] = RT_Table(stack(3).texture, options, lambda(i1),  strcat(stack(3).name, '_fi'));
        end
        
    elseif strcmp(stack(3).type, 'RCWA')
        disp('   Back surface RCWA')
        [C_R, C_T, C_A, angle_vector, C_summary] = RCWA(SC_3_1, grating_3_1, options, i1);
        
        %elseif stack(3).type == 'TMM'
        %    [C_R, C_T] = TMM(stack(2).layer, stack(3).layer, substrate, options, lambda(i1));
        %end
         elseif strcmp(stack(3).type, 'TMM')    
       disp('   Back surface TMM (planar surface)')
        [C_R, C_T, C_A, angle_vector, C_summary] = TMM(bulk, stack(3).layer, substrate, options, lambda(i1));
          
    end
    % saving matrices
    angle_vector = angle_vector(1:length(angle_vector(:,1))/2, 2:3);
    disp('   Done, saving back matrix')
    C = [angle_vector C_R; ...
        0 0 C_A];
    mat = C;
    parsave(strcat(options.name, '/C_', stack(3).name, num2str(lambda(i1)), 'nm.mat'), mat);
    %save(strcat(options.name, '/C_', num2str(lambda(i1)), 'nm.mat'), 'C');
    mat = C_summary;
    parsave(strcat(options.name, '/C_', stack(3).name, num2str(lambda(i1)), 'nm_sum.mat'), mat);
    %    save(strcat(options.name, '/C_', num2str(lambda(i1)), 'nm_sum.mat'), 'C_summary');
    
    T_back = C_T;
    else disp('Back matrix already calculated')
    end
    % D will be calculated by the OPTOS code 
end


end