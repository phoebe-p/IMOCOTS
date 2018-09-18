function [tr, inv_tr, options] = make_texture(texture_string, options)

if strcmp(texture_string, 'regular_pyramids_i')
    tr = regular_pyramids_i();
    if isfield(options, 'opening_angle')
        inv_tr = regular_pyramids_u(options.opening_angle);
    else
        inv_tr = regular_pyramids_u();
    end
    options.usemean = false;
    
elseif strcmp(texture_string, 'regular_pyramids_u')
    tr = regular_pyramids_u();
    if isfield(options, 'opening_angle')
        inv_tr = regular_pyramids_i(options.opening_angle);
    else
        inv_tr = regular_pyramids_i();
    end
    options.usemean = false;
    
elseif strcmp(texture_string, 'random_pyramids_i')
    tr = random_pyramids_i();
    inv_tr = random_pyramids_u();
    options.usemean = true;
    
elseif strcmp(texture_string, 'random_pyramids_u')
    tr = random_pyramids_u();
    inv_tr = random_pyramids_i();
    options.usemean = true;
    
end