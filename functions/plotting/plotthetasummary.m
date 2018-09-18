figure
title(strcat('Absorption (\lambda = ', num2str(wavelength), ' nm)'))
hold on
plot(sin_a_v, C_summary(end, 1:100), '-r', 'DisplayName', 'incident from front)')
legend('show')

figure;
title(strcat('\lambda = ', num2str(wavelength), ' nm'))
subplot(2, 1, 1)
pcolor(sin_a_v, sin_a_v, C_summary(1:100, 1:100)); shading flat;
h = gca;
set(h, 'YDir', 'reverse');
colorbar
xlabel('sin(\theta_{in})')
ylabel('sin(\theta_{out})')
title('Reflection')

subplot(2, 1, 2)
pcolor(sin_a_v, sin_a_v, C_summary(101:200, 1:100)); shading flat;
h = gca;  % Handle to currently active axes
h = gca;
set(h, 'YDir', 'reverse');
colorbar
xlabel('sin(\theta_{in})')
ylabel('sin(\theta_{out})')
title('Transmission') 
