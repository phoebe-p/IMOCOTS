function [stratum] = build_triangle_grid(layer_params, d, grid_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

h=layer_params.thickness;
strip_t = layer_params.gridline_width/layer_params.N_strips;
offset = 0.5*strip_t/tand(60);
widths = linspace(layer_params.gridline_width-2*offset, 2*offset, layer_params.N_strips);

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
    block.pmt_index=grid_pmt_index;
    stripe.block{1}=block;
    block.c2=1;
    block.pmt_index=cladding_pmt_index;
    stripe.block{2}=block;
    
    stratum.stripe{i1}=stripe;
end


% clear stripe
% stripe.type=1; % inhomogeneous
% stripe.c1=10/d;
% clear block
% block.c2=30/d;
% block.pmt_index=grid_pmt_index;
% stripe.block{1}=block;
% block.c2=1;
% block.pmt_index=cladding_pmt_index;
% stripe.block{2}=block;
%
% stratum.stripe{1}=stripe;
%
% clear stripe
% stripe.type=1; % inhomogeneous
% stripe.c1=20/d;
% clear block
% block.c2=20/d;
% block.pmt_index=grid_pmt_index;
% stripe.block{1}=block;
% block.c2=1;
% block.pmt_index=cladding_pmt_index;
% stripe.block{2}=block;
%
% stratum.stripe{2}=stripe;
%
% clear stripe
% stripe.type=1; % inhomogeneous
% stripe.c1=30/d;
% clear block
% block.c2=10/d;
% block.pmt_index=grid_pmt_index;
% stripe.block{1}=block;
% block.c2=1;
% block.pmt_index=cladding_pmt_index;
% stripe.block{2}=block;
%
% stratum.stripe{3}=stripe;
%
% clear stripe
% stripe.type=0;
% stripe.c1 = 1;
% stripe.pmt_index = cladding_pmt_index;
%
% stratum.stripe{4}=stripe;

clear stripe
stripe.type=0;
stripe.c1 = 1;
stripe.pmt_index = cladding_pmt_index;

stratum.stripe{end+1}=stripe;

end

