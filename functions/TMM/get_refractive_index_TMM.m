function [medium]=get_refractive_index_TMM(medium, wavelengths)
for j1=1:length(medium)
    if isempty(medium(j1).material_filename)
        if length(medium(j1).refractive_index) == 1
            medium(j1).refractive_index=medium(j1).refractive_index*ones(1,length(wavelengths));
        end
    else
        nk=load(medium(j1).material_filename);
        [~, unique_lambda] = unique(nk(:,1));
        nk = nk(unique_lambda, :);
        nk(:,1)=nk(:,1)*1000; %convert from um to nm
        medium(j1).refractive_index=interp1(nk(:,1),nk(:,2),wavelengths, 'linear', 'extrap') + ...
                1i*interp1(nk(:,1),nk(:,3),wavelengths, 'linear', 'extrap');
    end
end

end