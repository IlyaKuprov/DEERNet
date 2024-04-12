% Gunnar Jeschke's test database containing 75 DEER data sets with 
% known ground truth distance distributions. Sets 1-25 were designed
% to present challenging tests with respect to resolution of several
% Gaussian peaks, superposition of narrow and broad components, and
% shape recognition of broad asymmetric distributions. Sets 26-75 
% were derived from an ensemble model of the RNA-binding protein 
% PTBP1 in complex with an internal ribosome entry site of encephalo-
% myocarditis virus by selecting spin label site pairs from a set of
% labeling sites used in an ongoing experimental study. Files:
% 
% 1. Ground truth distance distributions
%    File names   : noise_test_pcr_?.dat (? = 1... 75)
%    First column : distance axis (Angstroem)
%    Second column: probability density, normalized to unity sum
% 
% 2. Low-noise simulated DEER data
%    File names   : low_noise_?.dat (? = 1... 75)
%    First column : time axis (nanoseconds)
%    Second column: real part of simulated DEER data
%    Third column : imaginary part of simulated DEER data (noise only)
% 
% 3. High-noise simulated DEER data
%    File names   : high_noise_?.dat (? = 1... 75)
%    First column : time axis (nanoseconds)
%    Second column: real part of simulated DEER data
%    Third column : imaginary part of simulated DEER data (noise only)
% 
% 4. Simulated DEER data with defined signal-to-noise ratios
%    File names   : SNR_$_set_?.dat ($ = 4,8,10,16,32,50,100; ? = 1... 75)
%    First column : time axis (nanoseconds)
%    Second column: real part of simulated DEER data
%    Third column : imaginary part of simulated DEER data (noise only)
%
% Any user of DEER data processing packages would do well to meditate
% on the graphical output of this script.
%
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function example_set_noise()

% Loop over data instances
for n=1:75
    
    % Process low-noise dataset
    raw_data=load(['data_noise/low_noise_' num2str(n) '.dat'],'-ASCII');
    [~,maxpt]=max(raw_data(:,2));
    deer_trace=raw_data(maxpt:end,2);
    time_axis=raw_data(maxpt:end,1);
    time_axis=time_axis-time_axis(1); 
    deernet(deer_trace,1e-9*time_axis);
    set(gcf,'Name',['Data set ' num2str(n) ', low noise'],...
            'NumberTitle','off'); drawnow();
    
    % Degrade the SNR
    for k=[100 50 32 16 10 8 4]
        
        % Process graded noise dataset
        raw_data=load(['data_noise/SNR_' num2str(k) '_set_' num2str(n) '.dat'],'-ASCII');
        deer_trace=raw_data(maxpt:end,2);
        time_axis=raw_data(maxpt:end,1);
        time_axis=time_axis-time_axis(1);
        deernet(deer_trace,1e-9*time_axis);
        set(gcf,'Name',['Data set ' num2str(n) ', SNR ' num2str(k)],...
                'NumberTitle','off'); drawnow();
        
    end
    
    % Process high-noise dataset
    raw_data=load(['data_noise/high_noise_' num2str(n) '.dat'],'-ASCII');
    deer_trace=raw_data(maxpt:end,2);
    time_axis=raw_data(maxpt:end,1);
    time_axis=time_axis-time_axis(1); 
    deernet(deer_trace,1e-6*time_axis);
    set(gcf,'Name',['Data set ' num2str(n) ', high noise'],...
            'NumberTitle','off'); drawnow();
    
    % Close figure windows
    pause(5); close('all');
    
end

end

