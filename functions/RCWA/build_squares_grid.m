function [stratum] = build_squares_grid(layer_params, d, grid_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

square = layer_params.gridline_width;
h=layer_params.thickness;

clear stratum
stratum.type=2; % biiperiodic stratum
stratum.thick=h; % stratum thickness
% The following h11 ... h22 spec defines the stratum's period vectors
% (GD-Calc.pdf, equation 3.18) and determines the stripe orientation.

% stratum vectors are parallel to grating vectors
stratum.h11=1;
stratum.h12=0;
stratum.h21=0;
stratum.h22=1;

clear stripe
stripe.type=1; % inhomogeneous
stripe.c1=square/d;
clear block
block.c2=square/d;
block.pmt_index=grid_pmt_index;
stripe.block{1}=block;
block.c2=1;
block.pmt_index=cladding_pmt_index;
stripe.block{2}=block;
stratum.stripe{1}=stripe;
clear stripe
stripe.type=0; % homogeneous
stripe.c1=1;
stripe.pmt_index=cladding_pmt_index;
stratum.stripe{2}=stripe;

end

