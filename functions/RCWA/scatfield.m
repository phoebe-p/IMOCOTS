function [grating, SC, R, T, inc_field] = scatfield(SC, grating, index)
wavelength=SC.wavelengths(index);
% Define the indicent field.
clear inc_field
inc_field.wavelength=wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1=-cos(theta)/wavelength.
sup_n = real(sqrt(cell2mat(grating.pmt(grating.pmt_sup_index))));
inc_field.f2=sup_n*sin(SC.theta)*cos(SC.phi)/wavelength;
inc_field.f3=sup_n*sin(SC.theta)*sin(SC.phi)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
order=[];
for m2=-SC.m_max:SC.m_max
    order(end+1).m2=m2;
    order(end).m1=-SC.m_max:SC.m_max;
end

% Run the diffraction calculations.
if any([SC.layer.type]==2)
    %save('g_before', 'grating')
    [param_size,scat_field,inc_field, grating]=gdc(grating,inc_field,order);
    %save('g_after', 'grating')
else
    [param_size,scat_field,inc_field]=gdc(grating,inc_field,order);
end
if isempty(scat_field)
    disp('Interrupted by user.');
    return
end

%save('scat_field','scat_field')
%save('inc_field', 'inc_field')
%save('grating','grating')

% Compute the diffraction efficiencies.
[R,T]=gdc_eff(scat_field,inc_field);
R_unpolarized=([R.eff1]+[R.eff2])/2;
T_unpolarized=([T.eff1]+[T.eff2])/2;
%store the total reflection, absorption and transmission for s an p
%incidence. The tranmsission is into the bottom cell. The absorption is
%everywhere else (middle and top cells, window, ARC and NPs)
SC.RAT.Rs(index)=sum([R.eff1]);
SC.RAT.Rp(index)=sum([R.eff2]);
SC.RAT.Rtotal(index)=sum(R_unpolarized);
SC.RAT.Ts(index)=sum([T.eff1]);
SC.RAT.Tp(index)=sum([T.eff2]);
SC.RAT.Ttotal(index)=sum(T_unpolarized);
SC.RAT.As(index)=1-sum([R.eff1]+[T.eff1]);
SC.RAT.Ap(index)=1-sum([R.eff2]+[T.eff2]);
SC.RAT.Atotal(index)=1-sum(R_unpolarized+T_unpolarized);


end
