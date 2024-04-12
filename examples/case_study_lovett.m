% DEER data from https://doi.org/10.1093/nar/gkaa086 with lots of
% tail artefacts that DEERNet diligently ignores.
%
% i.kuprov@soton.ac.uk

function case_study_lovett()

% Enumerate the data files
data_files=dir('data_lovett/*/*.DTA');

% Run the processing
for n=1:numel(data_files)
    
    % Load the data
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))]);
    
     % Run the processing
    deernet(deer_trace,time_axis,expt='deer'); 
    
    % Set the plot title
    set(gcf,'Name',data_files(n).name(1:(end-4)),...
            'NumberTitle','off'); drawnow();
    
end

end

