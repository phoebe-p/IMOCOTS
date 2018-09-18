function [opt_list] = plot_optical_constants(superstrate, stack, substrate, wavelengths)

opt_list(1).material.filename = superstrate;
filenames =  {superstrate};

for i1 = 1:length(stack)
    for j1 = 1:length(stack(i1).layer)
        if ~any(strcmp(filenames,stack(i1).layer(j1).mat{1}))
        opt_list(end+1).material.filename = stack(i1).layer(j1).mat{1};
        filenames(end+1) = stack(i1).layer(j1).mat;
        end
    end
end

opt_list(end+1).material.filename = substrate;

[opt_list]=get_refractive_index(opt_list, wavelengths);

figure
hold on
cols = lines(length(opt_list));
for i1 = 1:length(opt_list)
    yyaxis left
    plot(wavelengths/1000, real(opt_list(i1).material.pmt), '-', 'DisplayName', opt_list(i1).material.filename, 'Color', cols(i1,:), 'LineWidth', 1)
    yyaxis right
    h = plot(wavelengths/1000, imag(opt_list(i1).material.pmt), '--', 'Color', cols(i1,:), 'LineWidth', 1);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
end
set(gca, 'XScale', 'log');
yyaxis left
ylabel('n')
xlabel('Wavelength (um)')
l = legend('show');
set(l, 'Interpreter', 'none')
yyaxis right
set(gca, 'YScale', 'log');
ylabel('\kappa')
end

