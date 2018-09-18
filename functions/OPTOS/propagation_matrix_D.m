function [D_si] = propagation_matrix_D(cell_thickness,alpha, points)
disp('Calculation of propagation matrix...')
% calculate propagation matrix D
D_si_diag = zeros(size(points,2),1);
%size(D_si_diag)
for i1 = 1:size(points,2)
    theta = points(1,i1);
    D_si_diag(i1) = exp(-alpha.*cell_thickness./cos(theta*pi/180));
end
D_si_diag(size(points,2)-1:size(points,2)) = 1; % these entries represent the absorbed light
D_si = sparse(diag(D_si_diag,0));
disp('Calculation of propagation matrix done.')
end