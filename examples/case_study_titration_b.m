% Titration data from Nina Kubatova and Thomas Schmidt, 
% particulars in
%
%       https://doi.org/10.1073/pnas.2221036120
%
% nina.kubatova@nih.gov

function case_study_titration_b()

% Enumerate the data files
data_files=dir('data_titration\titration_cholate\*\3us.DTA');

% Run the processing
for n=1:numel(data_files)

    % Load the data
    [deer_trace,time_axis]=elexsys2deernet([data_files(n).folder filesep ...
                                            data_files(n).name(1:(end-4))]);

    % Run the processing 
    deernet(deer_trace,time_axis,expt='deer',bg_dim_range=[2.0,3.5]); 

    % Set the plot title
    set(gcf,'Name',data_files(n).folder,'NumberTitle','off'); drawnow();

end

end

