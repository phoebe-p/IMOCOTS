function [ ] = plot_grating( grating , d)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

clear pmt_display
for j1 = 1:grating.pmt_sup_index
    pmt_display(j1).name='';
    pmt_display(j1).color=[1,1,1];
    pmt_display(j1).alpha=1;
    grating.pmt{j1}=1;
end

pmt_display(3).color=[0.8 0.8 0.8]; % Ti
pmt_display(2).color=[1 0 0]; % Air
%pmt_display(4).color=[1, 0.8, 0]; %Au
%pmt_display(5).alpha=0; % GaAs?
%pmt_display(grating.pmt_sub_index).color=[0, 1, 0];

x_limit=[0,-1.5*d,-1.3*d;grating.total_thickness,1.5*d,1.3*d];
%x_limit=[0,0,0;grating.total_thickness,1*d,1*d];
gdc_plot(grating,1,pmt_display,x_limit);
end

