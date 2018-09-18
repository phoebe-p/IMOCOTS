function [v0,B,C,points] = load_matrices(wavelength, options, frontname, backname)
disp('Matrix loading...')

%front matrix (for light coming from bulk of SC)
matrix_name_B = strcat('results/', options.name, '/B_', frontname, num2str(wavelength), 'nm.mat');
B = load(matrix_name_B);
fn_B = fieldnames(B);
B = B.(fn_B{1});
% add a row to B to deal with back absorption:
B(:, end+1:end+2) = 0; % add columns of zeros (front & back absorption)

add_row = zeros(1, length(B(1,:)));
add_row(end) = 1;
B(end+1, :) = add_row;
B(end-1, end-1) = 1; % already absorbed light not affected
B = sparse(B);

%back  matrix
matrix_name_C =  strcat('results/', options.name, '/C_', backname, num2str(wavelength), 'nm.mat');
C = load(matrix_name_C);
fn_C = fieldnames(C);
C = C.(fn_C{1});
C(:, end+1:end+2) = 0;

% move last row (absorption) of C down one row, insert row in between to
% deal with front absorption
add_row = zeros(1, length(C(1,:)));
add_row(end-1) = 1;
C = [C(1:end-1, :); add_row; C(end,:)];
C(end, end) = 1;
C = sparse(C);



% need to add an extra column to B and C to deal with absorption (all
% zeros)
%B(:, end+1) = 0;
%C(:, end+1) = 0;


% delete first two columns with angle information
B(:,1:2) = [];
C(:,1:2) = [];

% calculate v0 and points
matrix_name_A = strcat('results/', options.name, '/A_', frontname, num2str(wavelength), 'nm.mat');
[v0, points] = calc_v_0(matrix_name_A, options);
v0 = sparse(v0);

points = [points; 0 0; 0 0]'; 
% two additional rows: (end-1) is front absorption, end is back absorption

%if TETM == 'TE' 																						% for TE-polarized incoming light
%				v0 = v0(:,3);
%		elseif TETM == 'TM'																						% for TM-polarized incoming light
%						v0 = v0(:,4);
%			end
disp('Matrix loading done.')
