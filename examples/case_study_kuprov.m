% The six examples used in the paper by Worswick, Spencer, 
% Jeschke, and Kuprov:
%
%         https://doi.org/10.1126/sciadv.aat5218
%
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function case_study_kuprov()

% Sample I
expt_data=load('data_kuprov/sample_I_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),expt='deer'); drawnow();

% Sample II
expt_data=load('data_kuprov/sample_II_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),expt='deer'); drawnow();

% Sample III
expt_data=load('data_kuprov/sample_III_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),expt='deer'); drawnow();

% Sample IV
expt_data=load('data_kuprov/sample_IV_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),expt='deer'); drawnow();

% Sample V
expt_data=load('data_kuprov/sample_V_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),expt='deer'); drawnow();

% Sample VI
expt_data=load('data_kuprov/sample_VI_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),expt='deer'); drawnow();

end

