function [grating, SC] = build_grating(SC)
[SC.layer]=get_refractive_index(SC.layer,SC.wavelengths);
%if any of the layers are a DBR then interpolate refractive indices and calculate thickness for given centre wavelength
for j1=1:length(SC.layer)
    if SC.layer(j1).type == 3
        nk_DBR1=load(SC.layer(j1).material(1).filename);
        nk_DBR2=load(SC.layer(j1).material(2).filename);
        n_DBR1=interp1(nk_DBR1(:,1)*1000,nk_DBR1(:,2), SC.layer(j1).centre_wavelength);
        n_DBR2=interp1(nk_DBR2(:,1)*1000,nk_DBR2(:,2), SC.layer(j1).centre_wavelength);
        SC.layer(j1).material(1).thickness=SC.layer(j1).centre_wavelength/n_DBR1/4;
        SC.layer(j1).material(2).thickness=SC.layer(j1).centre_wavelength/n_DBR2/4;
    end
end
[grating, SC]=build_generalized_structure(SC);
end