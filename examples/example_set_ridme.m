% RiDME data from Gunnar Jeschke's group at ETH and Janet Lovett's group
% at the University of St Andrews. The latter dataset is from
%
%              https://doi.org/10.1007/s00723-021-01321-6
%
% jake.keeley@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% jel20@st-andrews.ac.uk
% i.kuprov@soton.ac.uk

function example_set_ridme()

% Enumerate the data files
data_files=dir('./data_ridme/*.DTA');

% Run the processing
for n=1:numel(data_files)
    
    % Load RIDME data, truncating the gating spikes (first five points)
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))],[5 0]);
                                        
    % Run the processing
    deernet(deer_trace,time_axis,expt='ridme'); drawnow();
    
    % Set the plot title
    set(gcf,'Name',data_files(n).name(1:(end-4)),'NumberTitle','off');
    
end

end

