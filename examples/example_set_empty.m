% Datasets that do not contain any actual dipolar modulation
% signal, provided by Laura Galazzo and Enrica Bordignon. 
%
% laura.galazzo@ruhr-uni-bochum.de

function example_set_empty()

% Enumerate the data files
data_files=dir('./data_empty/*.dat');

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
    deernet(deer_trace,time_axis); drawnow();
    
end

end

