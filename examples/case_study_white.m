% Ring test set up by Olav Schiemann and Gunnar Jeschke for 
% the DEER data acquisition and processing white paper:
%
%         https://doi.org/10.1021/jacs.1c07371
%
% gunnar.jeschke@phys.chem.ethz.ch
% schiemann@pc.uni-bonn.de

function case_study_white()

% Enumerate the data files
data_files=dir('./data_white/*.DTA');

% Run the processing
for n=1:numel(data_files)
    
    % Load the data
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))]);
                                        
    % Cut off excessive baseline tails
    switch data_files(n).name(1:(end-9))
        
        case 'sample1'
            
            deer_trace(time_axis>3e-6)=[];
            time_axis(time_axis>3e-6)=[];
            
        case 'sample2'
            
            deer_trace(time_axis>4e-6)=[];
            time_axis(time_axis>4e-6)=[];
            
        case 'sample3'
            
            deer_trace(time_axis>6e-6)=[];
            time_axis(time_axis>6e-6)=[];
                
    end
    
    % Lab B appears to have a gating problem
    switch data_files(n).name(1:(end-4))
    
        case 'sample2_labB'
            
            deer_trace=deer_trace(15:end);
            time_axis=time_axis(15:end);
            time_axis=time_axis-time_axis(1);
            
        case 'sample3_labB'
            
            deer_trace=deer_trace(10:end);
            time_axis=time_axis(10:end);
            time_axis=time_axis-time_axis(1);
            
        case 'sample4_labB' 
            
            deer_trace=deer_trace(15:end);
            time_axis=time_axis(15:end);
            time_axis=time_axis-time_axis(1);
        
    end
                                        
    % Run the processing
    deernet(deer_trace,time_axis,bg_dim_range=[2.0 3.5]);
    
    % Set the plot title
    set(gcf,'Name',data_files(n).name(1:(end-4)),...
            'NumberTitle','off'); drawnow();
    
end

end

