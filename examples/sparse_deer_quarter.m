% DEER data from the example set of DeerAnalysis, kindly
% provided by Gunnar Jeschke. Undersampled ex post facto
% to 1/4 of the points to emulate sparse acquisition.
%
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function sparse_deer_quarter()

% Enumerate the data files
data_files=dir('./data_deeran/*.DTA');

% Run the processing
for n=1:numel(data_files)
    
    % Load the data
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))]);
    
    % Resample to 512 points
    resamp_axis=linspace(time_axis(1),time_axis(end),512)';
    resamp_trace=interp1(time_axis,deer_trace,resamp_axis,'pchip');
    
    % Build sparse sampling mask
    mask=false(512,1); 
    mask([1:127 512])=true;
    mask=mask([1 randperm(510)+1 512]);

    % Apply sparse sampling
    resamp_trace(~mask)=NaN;
    
    % Run DEERNet
    deernet(resamp_trace,resamp_axis); drawnow;
    
end

end

