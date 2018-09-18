function [stack] = make_effective_medium(stack, start_layer, end_layer, n_layers, thickness, wavelength)
% make for TMM, but maybe should also be available for RCWA?
% need two materials and thickness, and number of layers.
% Aspnes, Thin Solid Films, 1989

frac_start = linspace(1, 0, n_layers);
frac_end = 1 - frac_start;
endpoints(1).material_filename = start_layer;
endpoints(2).material_filename = end_layer;
[opt_consts]=get_refractive_index_TMM(endpoints, wavelength);
layer_thickness = thickness/n_layers;

for i1 = 1:n_layers
    stack(end+1).refractive_index=frac_start(i1)*opt_consts(1).refractive_index + ...
        frac_end(i1)*opt_consts(2).refractive_index;
    stack(end).thickness=layer_thickness;          %nm
    stack(end).material_filename = '';
end


end