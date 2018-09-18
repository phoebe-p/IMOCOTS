function [grating, SC, R_w, T_w, inc_field]=optical(grating, SC)
layer_indices = needed_layers(SC);
R_w = cell(length(SC.wavelengths), 1);
T_w = cell(length(SC.wavelengths), 1);
for j1=1:length(SC.wavelengths)
    %set pmts for this wavelength
    grating.pmt={};
    grating.pmt{1}=SC.layer(1).material.pmt(j1);
    for j2=2:length(SC.layer)
        for j3=1:length(SC.layer(j2).material)
            grating.pmt{end+1}=SC.layer(j2).material(j3).pmt(j1);
        end        
    end
    %calculate scattered fields - internal and external
    [grating, SC, R_w{j1}, T_w{j1}, inc_field] = scatfield(SC, grating, j1);
    [grating, SC] = Poynting_abs(SC, grating, j1, layer_indices);
end
end

function[grating, SC] = Poynting_abs(SC, grating, index, layer_indices)
        calc_layer = layer_indices.calc_layer;
        t = layer_indices.t;
        b = layer_indices.b;
        
        P_needed = sort(unique([b'; t']));
        P_data = table(calc_layer', t', b', 'VariableNames',{'i', 't', 'b'});
        for j = 1:length(P_needed)
            a = struct2table(grating.stratum{1, P_needed(j)}.full_field);
            ffE2s = table2array(a(:, 7));
            ffE3s = table2array(a(:, 8));
            ffH2s = table2array(a(:, 10));
            ffH3s = table2array(a(:, 11));
            
            ffE2p = table2array(a(:, 13));
            ffE3p = table2array(a(:, 14));
            ffH2p = table2array(a(:, 16));
            ffH3p = table2array(a(:, 17));

            P_needed(j, 2) = sum(real(ffE2s.*conj(ffH3s) - ffE3s.*conj(ffH2s)));
            P_needed(j, 3) = sum(real(ffE2p.*conj(ffH3p) - ffE3p.*conj(ffH2p)));
            % p polarization!
            
            %zeroth order:
            % zeroth order; 
            ind = find(table2array(a(:,1)) == 0 & table2array(a(:,2)) == 0);
            P_needed(j, 4) = real(ffE2s(ind).*conj(ffH3s(ind)) - ffE3s(ind).*conj(ffH2s(ind)));
            P_needed(j, 5) = real(ffE2s(ind).*conj(ffH3s(ind)) - ffE3s(ind).*conj(ffH2s(ind)));
            %%%%----%
            
        end
        findP_s = @(x) - P_needed(P_needed(:,1) == x, 2);
        findP_p = @(x) - P_needed(P_needed(:,1) == x, 3);
        findP_s_zeroth = @(x) - P_needed(P_needed(:,1) == x, 4);
        findP_p_zeroth = @(x) - P_needed(P_needed(:,1) == x, 5);
        %%

        P_data.P_abs_s = (arrayfun(findP_s, P_data.t) - arrayfun(findP_s, P_data.b))/cos(SC.theta);
        P_data.P_abs_p = (arrayfun(findP_p, P_data.t) - arrayfun(findP_p, P_data.b))/cos(SC.theta);
        P_data.P_abs_s_zeroth = (arrayfun(findP_s_zeroth, P_data.t) - arrayfun(findP_s_zeroth, P_data.b))/cos(SC.theta);
        P_data.P_abs_p_zeroth = (arrayfun(findP_p_zeroth, P_data.t) - arrayfun(findP_p_zeroth, P_data.b))/cos(SC.theta);
        for j2 = 1:length(P_data.i)
        SC.layer(P_data.i(j2)).P_abs_s(index) = P_data.P_abs_s(j2);
        SC.layer(P_data.i(j2)).P_abs_p(index) = P_data.P_abs_p(j2);
        SC.layer(P_data.i(j2)).P_abs(index) = (P_data.P_abs_s(j2) + P_data.P_abs_p(j2))/2;
        
        SC.layer(P_data.i(j2)).P_abs_s_zeroth(index) = P_data.P_abs_s_zeroth(j2);
        SC.layer(P_data.i(j2)).P_abs_p_zeroth(index) = P_data.P_abs_p_zeroth(j2);
        SC.layer(P_data.i(j2)).P_abs_zeroth(index) = (P_data.P_abs_s_zeroth(j2) + P_data.P_abs_p_zeroth(j2))/2;
        end
  
end




