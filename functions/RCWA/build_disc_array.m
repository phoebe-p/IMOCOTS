function [stratum] = build_disc_array(grating, layer_params , d, cylinder_pmt_index, cladding_pmt_index)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    N=layer_params.N_stripes;
    r=layer_params.radius;
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

    stratum.stripe=cell(1,2*N);
    % Define vertex coordinates for block-partitioned unit circle (see
    % GD-Calc.pdf, Figure 4). The j-th block vertex in the first quadrant
    % has coordinates [x(j),x(N+1-j)], j=1...N. x(j) is monotonic
    % decreasing with j.
    x=circle_partition(N);
    clear stripe block
    stripe.type=1; % inhomogeneous stripe
    % Construct the stratum. (Refer to Figures 3 and 4 in GD-Calc.pdf and
    % view the grating plot with a small N value to follow the construction
    % logic.)
    % The x2, x3 coordinate origin is at the center of a pillar.
    for n=1:N
        if n<N
            % The next stripe intercepts a row of pillars between the x3
            % coordinate limits -x(n)*r and -x(n+1)*r.
            stripe.c1=-x(n+1)*r/grating.d31;
        else
            % The N-th stripe is centered on the pillar axes, and its x3
            % limits are -x(N)*r and +x(N)*r.
            stripe.c1=x(N)*r/grating.d31;
        end
        % The first block defines the open space between adjacent pillars.
        % The block's x2 coordinate limits are x(N+1-n)*r-d and
        % -x(N+1-n)*r.
        block.c2=(-x(N+1-n)*r-stripe.c1*grating.d21)/d;
        block.pmt_index=cladding_pmt_index;
        stripe.block{1}=block;
        % The second block traverses the pillar interior. Its x2 coordinate
        % limits are -x(N+1-n)*r and +x(N+1-n)*r.
        block.c2=(x(N+1-n)*r-stripe.c1*grating.d21)/d;
        block.pmt_index=cylinder_pmt_index;
        stripe.block{2}=block;
        stratum.stripe{n}=stripe;
    end
    for n=2:N
        % The next stripe intercepts a row of pillars between the x3
        % coordinate limits x(N+2-n)*r and x(N+1-n)*r.
        stripe.c1=x(N+1-n)*r/grating.d31;
        % The first block defines the open space between adjacent pillars.
        % The block's x2 coordinate limits are x(n)*r-d and -x(n)*r.
        block.c2=(-x(n)*r-stripe.c1*grating.d21)/d;
        block.pmt_index=cladding_pmt_index;
        stripe.block{1}=block;
        % The second block traverses the pillar interior. Its x2 coordinate
        % limits are -x(n)*r and +x(n)*r.
        block.c2=(x(n)*r-stripe.c1*grating.d21)/d;
        block.pmt_index=cylinder_pmt_index;
        stripe.block{2}=block;
        stratum.stripe{N+n-1}=stripe;
    end
    clear stripe
    % The next stripe defines the open space between adjacent rows of
    % pillars. Its x3 coordinate limits are x(1)*r and grating.d31-x(1)*r.
    stripe.type=0; % homogeneous stripe
    stripe.c1=1-x(1)*r/grating.d31;
    stripe.pmt_index=cladding_pmt_index;
    stratum.stripe{2*N}=stripe;
    clear stripe
  

end

