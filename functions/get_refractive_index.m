function [medium]=get_refractive_index(medium, wavelengths)
for j1=1:length(medium)
    for j2=1:length(medium(j1).material)
        %  if length(medium(j1).material(j2).refractive_index)~=length(wavelengths)
        if isempty(medium(j1).material(j2).filename)
            medium(j1).material(j2).refractive_index = ...
                medium(j1).material(j2).refractive_index*ones(1,length(wavelengths));
        else
            nk=load(medium(j1).material(j2).filename);
            [~, unique_lambda] = unique(nk(:,1));
            nk = nk(unique_lambda, :);
            nk_interp=interp1(nk(:,1)*1000,nk(:,2),wavelengths, 'linear', 'extrap')+ ...
                1i*interp1(nk(:,1)*1000,nk(:,3),wavelengths, 'linear', 'extrap');
            medium(j1).material(j2).refractive_index=nk_interp;
        end
        medium(j1).material(j2).pmt = ...
            medium(j1).material(j2).refractive_index.*medium(j1).material(j2).refractive_index;
        %  end
    end
end

end