function [grating, SC] = make_grating(superstrate, stack_entry, substrate, options)
% define structure

% TYPES:
% 1: homogeneous
% 2: homogeneous with full-field output
% 3: Distributed Bragg Reflector (DBR)
% 4: Disk array
% 5: Biperiodic crossed grid ('Massiot')
% 6: Uniperiodic line grid
% 7: Array of square pillars (biperiodic)
% 8: Array of spheres
% 9: coordinate break (for e.g. tapered pillars)
% 10: OPTOS-paper style diagonal array of squares
% 11: triangles in triangular array
% 12: hexagons in triangular array

SC.wavelengths=options.wavelengths;

SC.period=stack_entry.grid_period;
%SC.N=80; % number of partition blocks per semicircle (for pillars)
SC.m_max=stack_entry.m_max; % maximum diffraction order index

%define layer stack--------------------------------------------------------
%build grating from bottom to top

%layer(1) is semi-infinite transmission medium
SC.layer(1).material.filename=substrate;
SC.layer(1).material.refractive_index=[];

% have not implemented DBR or disk array

for i1 = flip(1:length(stack_entry.layer)) % GD-Calc builds structures bottom to top!
    if strcmp(stack_entry.layer(i1).grid_type, 'homogeneous')
        SC.layer(end+1).type = 1;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'homogeneous_ff')
        SC.layer(end+1).type = 2;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'DBR')
        SC.layer(end+1).type = 3;
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'disk_array')
        SC.layer(end+1).type = 4;
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'crossed_grid')
        SC.layer(end+1).type = 5;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).gridline_width=stack_entry.layer(i1).line_width;    %nm
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'line_grid')
        SC.layer(end+1).type = 6;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).gridline_width=stack_entry.layer(i1).line_width;    %nm
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'square_array')
        SC.layer(end+1).type = 7;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).gridline_width=stack_entry.layer(i1).line_width;    %nm
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'sphere_array')
        SC.layer(end+1).type = 8;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).N=stack_entry.layer(i1).N;    %nm
        SC.layer(end).zeta =stack_entry.layer(i1).zeta; % in radians
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'coordinate_break')
        SC.layer(end+1).type = 9;
        SC.layer(end).shift_x2 = stack_entry.layer(i1).shift_x2; % half-period x2-shift
        SC.layer(j1).shift_x3 = stack_entry.layer(i1).shift_x23;
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'diag_square_array')
        SC.layer(end+1).type = 10;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).gridline_width=stack_entry.layer(i1).line_width;
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'diag_triangle_array')
        SC.layer(end+1).type = 11;
        SC.layer(end).N_strips = stack_entry.layer(i1).N_strips;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).gridline_width=stack_entry.layer(i1).line_width;
        
    elseif strcmp(stack_entry.layer(i1).grid_type, 'diag_hexagon_array')
        SC.layer(end+1).type = 12;
        SC.layer(end).N_strips = stack_entry.layer(i1).N_strips;
        SC.layer(end).material(1).filename=stack_entry.layer(i1).mat{1};      %grid material
        SC.layer(end).material(1).refractive_index=[];
        SC.layer(end).material(2).filename=stack_entry.layer(i1).mat{2};       %cladding material
        SC.layer(end).material(2).refractive_index=[];
        SC.layer(end).thickness=stack_entry.layer(i1).thickness; %nm
        SC.layer(end).gridline_width=stack_entry.layer(i1).line_width;
    else
        error('invalid grating type for RCWA')
    end
    
end
SC.layer(end+1).material.filename=superstrate;
SC.layer(end).material.refractive_index=[];


[grating, SC] = build_grating(SC);
end

function [grating, SC] = build_grating(SC)
[SC.layer]=get_refractive_index(SC.layer,SC.wavelengths);
%if any of the layers are a DBR then interpolate refractive indices and calculate thickness for given centre wavelength
for j1=1:length(SC.layer)
    if SC.layer(j1).type == 3
        nk_DBR1=load(SC.layer(j1).material(1).filename);
        nk_DBR2=load(SC.layer(j1).material(2).filename);
        n_DBR1=interp1(nk_DBR1(:,1)*1000,nk_DBR1(:,2), SC.layer(j1).centre_wavelength);
        n_DBR2=interp1(nk_DBR2(:,1)*1000,nk_DBR2(:,2), SC.layer(j1).centre_wavelength);
        SC.layer(j1).material(1).thickness=SC.layer(j1).centre_wavelength/n_DBR1/4;
        SC.layer(j1).material(2).thickness=SC.layer(j1).centre_wavelength/n_DBR2/4;
    end
end
[grating, SC]=build_generalized_structure(SC);
end
