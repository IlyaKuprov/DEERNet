% Very broad distance distribution datasets provided by
% Gunnar Jeschke's lab, data from Irina Ritch. Figure 4
% from:
%
%       https://doi.org/10.3389/fmolb.2021.636599
%
% gunnar.jeschke@phys.chem.ethz.ch

function example_set_broad()

% Enumerate the data files
data_files=dir('./data_broad/*.dat');

% Run the processing
for n=1:numel(data_files)
    
    % Load each file
    raw_data=load([data_files(n).folder filesep data_files(n).name]); 
    time_axis=1e-6*raw_data(:,1); deer_trace=raw_data(:,2);
    
    % Chop off negative times
    deer_trace=deer_trace(time_axis>=0); 
    time_axis=time_axis(time_axis>=0);
    time_axis=time_axis-time_axis(1);
    
    % Run the processing
    deernet(deer_trace,time_axis); 
    
    % Set the plot title
    set(gcf,'Name',data_files(n).name(1:(end-4)),...
            'NumberTitle','off'); drawnow();
    
end

end

