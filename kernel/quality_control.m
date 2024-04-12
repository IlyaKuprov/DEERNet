% DEERNet quality control - heuristic tests designed to prevent
% unskilled users from shooting themselves in the foot. Syntax:
%
%               dataset=quality_control(dataset)
%
% This is an internal function, direct calls are discouraged.
%
% i.kuprov@soton.ac.uk 
%
% #NWIKI

function dataset=quality_control(dataset) 

% Check consistency
disp(newline); grumble(dataset);

% Loop over input traces
for n=1:dataset.ntraces
    
    % Do all nets return positive modulation depths?
    all_nets_positive_md=all(dataset.mdpths(1,n,:)>0,'all');
    
    % Is the average modulation depth sensible?
    decent_average_md=(mean(dataset.mdpths(1,n,:),'all')>0.005);

    % Is the uncertainty sensible?
    signal_norm=norm(dataset.dist_av(:,n),2);
    spread_norm=norm(dataset.dist_ub(:,n)-dataset.dist_av(:,n),2);
    decent_uncertainty=(spread_norm/signal_norm<0.6);
        
    % Check for nonsensical situations
    if (~all_nets_positive_md)||(~decent_average_md)||(~decent_uncertainty)
        dataset.dist_av(:,n)=NaN; 
        dataset.dist_lb(:,n)=NaN; % Prevent distance distributions
        dataset.dist_ub(:,n)=NaN; % from being returned
    end
    
    % Do we have a peak at the long edge?
    long_edge_peak=(dataset.dist_ub(end,n)>0.20*max(dataset.dist_av(:,n)));
    
    % Warn the user
    if long_edge_peak
        fprintf(['WARNING: significant distance density at the long edge;\n' ...
                 '         you are pushing your luck, consider recording\n'  ...
                 '         a longer trace and cutting off end artefacts.\n\n']);
    end
    
    % Do we have an unstable background estimate?
    bg_head_wobble=(dataset.backgs_ub(1,n)-dataset.backgs_lb(1,n))/dataset.backgs_av(1,n);
    bg_tail_wobble=(dataset.backgs_ub(end,n)-dataset.backgs_lb(end,n))/dataset.backgs_av(end,n);
    wobbly_background=(bg_head_wobble>0.40*dataset.mdpths_av(1,n))|...
                      (bg_tail_wobble>0.20*dataset.mdpths_av(1,n));
    
    % Warn the user
    if wobbly_background
        fprintf(['WARNING: networks disagree on the background signal;\n' ...
                 '         modulation depth cannot be estimated.\n\n']);
    end
    
    % Do not return unreliable backgrounds and modulation depths
    if long_edge_peak||wobbly_background||(~all_nets_positive_md)||(~decent_average_md)
        dataset.backgs_av(:,n)=NaN; dataset.backgs_lb(:,n)=NaN;
        dataset.backgs_ub(:,n)=NaN; dataset.mdpths_av(1,n)=NaN;
        dataset.mdpths_st(1,n)=NaN;
    end
    
end

% Clean up the output
dataset=rmfield(dataset,{'backgs','retros','mdpths'});

end

% Consistency enforcement
function grumble(dataset)
if ~isfield(dataset,{'ntraces','mdpths','dist_av','dist_lb',...
                     'dist_ub','backgs_ub','backgs_lb','backgs_av',...
                     'mdpths_av','mdpths_st','backgs','retros','mdpths'})
    error('the dataset object appears to be corrupted.');
end
end

% Dr George W. Beadle
% The University of Chicago
% Office of the President
% Chicago, Illinois
%
% Dear George,
%
% Yours is the first honorary degree that I have been offered, 
% and I thank you for considering me for such an honor.
% 
% However, I remember the work I did to get a real degree at 
% Princeton and the guys on the same platform receiving honorary
% degrees without work - and felt an "honorary degree" was a de-
% basement of the idea of a "degree which confirms certain work
% has been accomplished". It is like giving an "honorary electri-
% cians license". I swore then that if by chance I was ever offe-
% red one I would not accept it.
% 
% Now at last (twenty-five years later) you have given me a chan-
% ce to carry out my vow.
% 
% So thank you, but I do not wish to accept the honorary degree
% you offered.
% 
% Sincerely yours,
% Richard P. Feynman

