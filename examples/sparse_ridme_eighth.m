% RiDME data from Gunnar Jeschke's group at ETH and Janet 
% Lovett's group at the University of St Andrews. The lat-
% ter dataset is from
%
%        https://doi.org/10.1007/s00723-021-01321-6
%
% Undersampled ex post facto to 1/8 of the points to emu-
% late sparse acquisition.
%
% jake.keeley@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% jel20@st-andrews.ac.uk
% i.kuprov@soton.ac.uk

function sparse_ridme_eighth()

% Enumerate the data files
data_files=dir('./data_ridme/*.DTA');

% Run the processing
for n=1:numel(data_files)
    
    % Load RIDME data, truncating the gating spikes (first five points)
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))],[5 0]);

    % Resample to 512 points
    resamp_axis=linspace(time_axis(1),time_axis(end),512)';
    resamp_trace=interp1(time_axis,deer_trace,resamp_axis,'pchip');
    
    % Build sparse sampling mask
    mask=false(512,1); 
    mask([1:63 512])=true;
    mask=mask([1 randperm(510)+1 512]);

    % Apply sparse sampling
    resamp_trace(~mask)=NaN;
                                        
    % Run the processing
    deernet(resamp_trace,resamp_axis,expt='ridme'); drawnow();
    
    % Set the plot window title to the file name
    set(gcf,'Name',data_files(n).name(1:(end-4)),'NumberTitle','off');
    
end

end

