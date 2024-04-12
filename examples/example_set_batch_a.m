% Batch processing of titration data with a common baseline. The 
% data was acquired by Nina Kubatova and Thomas Schmidt:
%
%            https://doi.org/10.1073/pnas.2221036120
%
% nina.kubatova@nih.gov

function example_set_batch_a()

% Enumerate the data files
data_files=dir('data_titration\titration_chs\*\3us.DTA');

% Load the data
for n=1:numel(data_files)
    [deer_trace{n},time_axis{n}]=elexsys2deernet([data_files(n).folder filesep ...
                                                  data_files(n).name(1:(end-4))]); %#ok<AGROW>
end

% Truncate the data to line up the time axis
[shortest_set,idx]=min(cellfun(@numel,time_axis));
for n=1:numel(deer_trace)
    deer_trace{n}=deer_trace{n}(1:shortest_set);
end
time_axis=time_axis{idx}; deer_trace=cell2mat(deer_trace);

% Run the processing 
deernet(deer_trace,time_axis,expt='deer',bg_dim_range=[2.0,3.5]); 

end

