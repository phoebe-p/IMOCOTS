% set default options, add relevant folders to path
addpath('nk')
addpath(genpath('functions'))
warning('off')

options = [];
wavelengths = linspace(log(350), log(16*1000), 20);
wavelengths = floor(exp(wavelengths));
options.angle_in = 0;
options.wavelengths = wavelengths;
options.n_angle_bins = 100;
options.theta_max = 90;
options.phi_max = 90;
options.c_azimuth = 0.25;
options.n_scans = 5;
options.n_rays = 50000;
options.mirror = true;
options.theta_in = 0;
options.phi_in = 'uniform';
options.parallel = true;
options.incident_absorption = true;
options.power_threshold = 1e-18;
options.lookuptable_step = 0.1;

iter_or_geom = 'iter';	