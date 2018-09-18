function [grating, SC]=build_generalized_structure(SC)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

clear grating

grating.d21=0;

% if isfield(SC, 'period') == 1
%     grating.d31=SC.period;
%     grating.d22=SC.period;
% else
%     grating.d31=1;
%     grating.d22=1;
% end

grating.d31=SC.period;
grating.d22=SC.period;

grating.d32=0;

material_index=1;
layer_index=0;
total_thickness=0;
for j1=2:length(SC.layer)
    if SC.layer(j1).type == 1  %homogenous stratum
        material_index=material_index+1;
        layer_index=layer_index+1;
        clear stratum
        stratum.type=0;
        stratum.thick=SC.layer(j1).thickness;
        stratum.pmt_index=material_index;
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 2   %homogenous stratum with full-field calculation
        material_index=material_index+1;
        SC.layer(j1).grating_layer_indices=[layer_index+1, layer_index+2];
        %first layer has zero thickness
        layer_index=layer_index+1;
        clear stratum
        stratum.type=0;
        stratum.thick=0;
        stratum.pmt_index=material_index;
        stratum.full_field = [];
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        %second layer has nonzero thickness
        layer_index=layer_index+1;
        clear stratum
        stratum.type=0;
        stratum.thick=SC.layer(j1).thickness;
        stratum.pmt_index=material_index;
        stratum.full_field = [];
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 3  %DBR
        material_index=material_index+2;
        for j2=1:SC.layer(j1).n_layers
            layer_index=layer_index+1;
            clear stratum
            stratum.type=0; % homogenous stratum
            if mod(j2,2)==1
                stratum.thick=SC.layer(j1).material(1).thickness;
                stratum.pmt_index=material_index-1;
                total_thickness=total_thickness+stratum.thick;
            else
                stratum.thick=SC.layer(j1).material(2).thickness;
                stratum.pmt_index=material_index;
                total_thickness=total_thickness+stratum.thick;
            end
            grating.stratum{layer_index}=stratum;
        end
        SC.layer(j1).n_strata = SC.layer(j1).n_layers;
    elseif SC.layer(j1).type == 4  %periodic array of disks
        material_index=material_index+2;
        layer_index=layer_index+1;
        clear stratum
        [stratum] = build_disc_array(grating, SC.layer(j1) ,  SC.period, material_index-1,material_index);
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 5   %crossed grid (massiot)
        material_index=material_index+2;
        layer_index=layer_index+1;
        clear stratum
        [stratum] = build_crossed_grid(SC.layer(j1), SC.period, material_index-1,material_index);
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 6 % line grid
        material_index = material_index + 2;
        layer_index = layer_index + 1;
        clear stratum
        [stratum] =  build_line_grid(SC.layer(j1), SC.period, material_index-1,material_index);
        grating.stratum{layer_index} =  stratum;
        total_thickness = total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 7   %square grid
        material_index=material_index+2;
        layer_index=layer_index+1;
        clear stratum
        [stratum] = build_squares_grid(SC.layer(j1), SC.period, material_index-1,material_index);
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 8 % spherical nanoparticles
        material_index=material_index+2;
        %grating = build_spheres(SC.layer(j1), SC.period, material_index-1, material_index, grating);
        
        grid_pmt_index = material_index - 1;
        cladding_pmt_index = material_index;
        alt_stripe=1; % selector for alternate stripe orientation (1, 2, or 3)
        N = SC.layer(j1).N; % number of rectangular strips per semicircle at full sphere circumference
        radius = SC.layer(j1).thickness;
        zeta = SC.layer(j1).zeta; % angle between lattice vectors. zeta = 0 is a square grid of spheres
        d = SC.period;
        
        %d=400; % grating period
        
        x=circle_partition(N);
        
        % Construct grating.
        % Define the x2 and x3 projections of the first grating period
        % (d21,d31) and second grating period (d22,d32). The second period is
        % parallel to the x2 axis, and the first period is at an angle zeta to
        % the x3 axis.
        
        % THIS WILL OVERRIDE ANYTHING SET PREVIOUSLY
        grating.d21=d*sin(zeta);
        grating.d31=d*cos(zeta);
        grating.d22=d;
        grating.d32=0;
        
        heights = radius*abs(diff([x 0]));
        
        radii = radius*x;
        radii = sort(radii);
        radii = [radii fliplr(radii(1:(end-1)))];
        heights = [heights(1:(end-1)) 2*heights(end) fliplr(heights(1:(end-1)))];
        
        Ns = [1:(N-1) N (N-1):-1:1];
        
        for i1 = 1:length(heights)
            N = Ns(i1);
            x = circle_partition(N);
            r = radii(i1);
            h = heights(i1); % grating height
            layer_index = layer_index + 1;
            clear stratum
            stratum.type=2; % biperiodic stratum
            stratum.thick=h; % stratum thickness
            % The following h11 ... h22 spec defines the stratum's period vectors
            % (GD-Calc.pdf, equation 3.18) and determines the stripe orientation.
            if alt_stripe==1
                % Stripes are parallel to [grating.d22,grating.d32].
                stratum.h11=1;
                stratum.h12=0;
                stratum.h21=0;
                stratum.h22=1;
            elseif alt_stripe==2
                % Stripes are parallel to [grating.d21,grating.d31].
                stratum.h11=0;
                stratum.h12=1;
                stratum.h21=1;
                stratum.h22=0;
            else % alt_stripe==3
                % Stripes are parallel to
                % [grating.d21,grating.d31]-[grating.d22,grating.d32].
                stratum.h11=1;
                stratum.h12=1;
                stratum.h21=0;
                stratum.h22=-1;
            end
            stratum.stripe=cell(1,2*N);
            % Define vertex coordinates for block-partitioned unit circle (see
            % GD-Calc.pdf, Figure 4). The j-th block vertex in the first quadrant
            % has coordinates [x(j),x(N+1-j)], j=1...N. x(j) is monotonic
            % decreasing with j.
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
                block.pmt_index=grid_pmt_index;
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
                block.pmt_index=grid_pmt_index;
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
            grating.stratum{layer_index}=stratum;
            total_thickness=total_thickness+stratum.thick;
            clear stratum
        end
        
        SC.layer(j1).n_strata = 2*SC.layer(j1).N - 1;
        
    elseif SC.layer(j1).type == 9 % coordinate break
        layer_index = layer_index + 1;
        clear stratum
        stratum.type=3; % coordinate break
        stratum.dx2=SC.layer(j1).shift_x2; % half-period x2-shift
        stratum.dx3=SC.layer(j1).shift_x3; % half-period x3-shift
        grating.stratum{layer_index} = stratum;
    elseif SC.layer(j1).type == 10   %diagonal square grid
        material_index=material_index+2;
        layer_index=layer_index+1;
        clear stratum
        [stratum] = build_diag_squares_grid(SC.layer(j1), SC.period, material_index-1,material_index);
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    elseif SC.layer(j1).type == 11   %diagonal triangle grid
        grating.d21=SC.period;
        grating.d31=0;
        grating.d22=cosd(60)*SC.period;
        grating.d32=sind(60)*SC.period;

        material_index=material_index+2;
        layer_index=layer_index+1;
        clear stratum
        [stratum] = build_triangle_grid(SC.layer(j1), SC.period, material_index-1,material_index);
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
        elseif SC.layer(j1).type == 12   %diagonal hexagon grid
        grating.d21=SC.period;
        grating.d31=0;
        grating.d22=cosd(60)*SC.period;
        grating.d32=sind(60)*SC.period;

        material_index=material_index+2;
        layer_index=layer_index+1;
        clear stratum
        [stratum] = build_hexagon_grid(SC.layer(j1), SC.period, material_index-1,material_index);
        grating.stratum{layer_index}=stratum;
        total_thickness=total_thickness+stratum.thick;
        SC.layer(j1).n_strata = 1;
    end
end

grating.pmt_sub_index=1; % substrate permittivity index
grating.pmt_sup_index=material_index+1; % superstrate permittivity index
grating.total_thickness=total_thickness;

end

