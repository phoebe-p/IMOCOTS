function [stratum] = build_diag_squares_grid(layer_params, d, grid_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

square = layer_params.gridline_width;
gap = (d-2*layer_params.gridline_width)/2;
h=layer_params.thickness;

clear stratum
stratum.type=2; % biiperiodic stratum
stratum.thick=h; % stratum thickness
% The following h11 ... h22 spec defines the stratum's period vectors
% (GD-Calc.pdf, equation 3.18) and determines the stripe orientation.

% stratum vectors are at 45degrees to grating vectors
stratum.h11=1;
stratum.h12=1;
stratum.h21=-1;
stratum.h22=1;

clear stripe
stripe.type=1; % inhomogeneous
stripe.c1=square/d;
clear block
block.c2=(square+gap)/d;
block.pmt_index=cladding_pmt_index;
stripe.block{1}=block;
block.c2=(2*square+gap)/d;
block.pmt_index=grid_pmt_index;
stripe.block{2}=block;
block.c2 = 1;
block.pmt_index=cladding_pmt_index;
stripe.block{3}=block;

stratum.stripe{1}=stripe;

clear stripe
stripe.type=0; % homogeneous
stripe.c1=(square+gap)/d;
stripe.pmt_index=cladding_pmt_index;
stratum.stripe{2}=stripe;

clear stripe
stripe.type=1; % inhomogeneous
stripe.c1=(2*square+gap)/d;
clear block
block.c2=square/d;
block.pmt_index=grid_pmt_index;
stripe.block{1}=block;
block.c2=1;
block.pmt_index=cladding_pmt_index;
stripe.block{2}=block;

stratum.stripe{3}=stripe;

clear stripe
stripe.type=0; % homogeneous
stripe.c1=1;
stripe.pmt_index=cladding_pmt_index;
stratum.stripe{4}=stripe;
end

