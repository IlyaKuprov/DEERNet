% DEER data from the example set of DeerAnalysis, kindly
% provided by Gunnar Jeschke.
%
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function example_set_deeran()

% Enumerate the data files
data_files=dir('./data_deeran/*.DTA');

% Run the processing
for n=1:numel(data_files)
    
    % Load the data
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))]);
    % Run the processing
    deernet(deer_trace,time_axis,expt='deer'); drawnow();
    
    % Set the plot title
    set(gcf,'Name',data_files(n).name(1:(end-4)),'NumberTitle','off');
    
end

end

