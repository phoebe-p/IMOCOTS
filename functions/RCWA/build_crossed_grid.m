function [stratum] = build_crossed_grid(layer_params, d, grid_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    w=layer_params.gridline_width;
    h=layer_params.thickness;
    
    clear stratum
    stratum.type=2; % biperiodic stratum
    stratum.thick=h; % stratum thickness
    % The following h11 ... h22 spec defines the stratum's period vectors
    % (GD-Calc.pdf, equation 3.18) and determines the stripe orientation.

    % Stripes are parallel to [grating.d22,grating.d32].
    stratum.h11=1;
    stratum.h12=0;
    stratum.h21=0;
    stratum.h22=1;

    stratum.stripe=cell(1,2);
    %first stripe is inhomogeneous
    clear stripe block
    stripe.type=1; % inhomogeneous stripe
    stripe.c1=-w/2/d;
    block.c2=-w/2/d;
    block.pmt_index=cladding_pmt_index;
    stripe.block{1}=block;
    block.c2=w/2/d;
    block.pmt_index=grid_pmt_index;
    stripe.block{2}=block;
    stratum.stripe{1}=stripe;
    
    %second stripe is homogeneous
    clear stripe
    stripe.type=0; % inhomogeneous stripe
    stripe.c1=w/2/d;
    stripe.pmt_index=grid_pmt_index;
    stratum.stripe{2}=stripe;
    clear stripe
    
  

end

