function [stratum] = build_hexagon_grid(layer_params, d, grid_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

h=layer_params.thickness;
strip_t = layer_params.gridline_width/layer_params.N_strips;
offset = strip_t*cosd(60);
widths = linspace(layer_params.gridline_width-offset, offset, layer_params.N_strips);

clear stratum
stratum.type=2; % biperiodic stratum
stratum.thick=h; % stratum thickness
% The following h11 ... h22 spec defines the stratum's period vectors
% (GD-Calc.pdf, equation 3.18) and determines the stripe orientation.

% stratum vectors are at 45degrees to grating vectors
stratum.h11=1;
stratum.h12=0;
stratum.h21=0;
stratum.h22=1;

for i1 = 1:layer_params.N_strips
    clear stripe
    stripe.type=1; % inhomogeneous
    stripe.c1=i1*strip_t/d;
    clear block
    block.c2=widths(i1)/d;
    block.pmt_index=cladding_pmt_index;
    stripe.block{1}=block;
  
    block.c2=(2*layer_params.gridline_width-offset)/d;
    block.pmt_index=grid_pmt_index;
    stripe.block{2}=block;
    
    stratum.stripe{i1}=stripe;
end

widths = linspace(2*layer_params.gridline_width-2*offset, layer_params.gridline_width, layer_params.N_strips);

for i1 = 1:layer_params.N_strips
    clear stripe
    stripe.type=1; % inhomogeneous
    stripe.c1=(layer_params.gridline_width+i1*strip_t)/d;
    clear block
    block.c2=widths(i1)/d;
    block.pmt_index=grid_pmt_index;
    stripe.block{1}=block;
    %block.c2=0.1;
    %block.pmt_index=cladding_pmt_index;
    %stripe.block{1}=block;
    block.c2=1;
    block.pmt_index=cladding_pmt_index;
    stripe.block{2}=block;
    
    stratum.stripe{layer_params.N_strips+i1}=stripe;
end


clear stripe
stripe.type=0;
stripe.c1 = 1;
stripe.pmt_index = cladding_pmt_index;

stratum.stripe{end+1}=stripe;

end

