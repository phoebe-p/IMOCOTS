function [stratum] = build_line_grid(layer_params, d, grid_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    w=layer_params.gridline_width;
    h=layer_params.thickness;
    
    clear stratum
    stratum.type=1; % uniperiodic stratum
    stratum.thick=h; % stratum thickness
    % The following h11 ... h22 spec defines the stratum's period vectors
    % (GD-Calc.pdf, equation 3.18) and determines the stripe orientation.

    % Stripes are parallel to [grating.d22,grating.d32].
    stratum.h11=1;
    stratum.h12=0;

    clear stripe
    stripe.c1=-w/2/d;
    stripe.pmt_index=cladding_pmt_index;
    stratum.stripe{1}=stripe;
    clear stripe
    stripe.c1=w/2/d;
    stripe.pmt_index=grid_pmt_index;
    stratum.stripe{2}=stripe;

end

