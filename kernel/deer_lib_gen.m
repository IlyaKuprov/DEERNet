% Generates a library of distance distributions and corresponding DEER
% or RIDME traces for use in neural network training. Syntax:
%
%               library=deer_lib_gen(file_name,parameters)
%
% Required fields of the parameters.* structure:
%
%    max_time          - DEER trace duration, seconds
% 
%    max_exch          - maximum exchange coupling, fraction of
%                        the maximum frequency representable on
%                        the current time discretisation grid
%
%    max_exch          - minimum exchange coupling, fraction of
%                        the maximum frequency representable on
%                        the current time discretisation grid
%
%    ntraces           - number of traces you wish to generate
% 
%    ndistmax          - maximum number of skewed gaussians 
%                        in the distance distribution
%
%    np_time           - number of digitisation points in the
%                        DEER trace
%
%    np_dist           - number of digitisation points in the
%                        distance distribution
%
%    npt_acq           - number of digitization points actually
%                        acquired in a sparsely sampled dataset;
%                        points are distributed randomly with
%                        uniform sampling probability
%
%    noise_lvl         - maximum RMS noise level as a fraction
%                        of the modulation depth (min is zero)
%
%    min_fwhm          - minimum FWHM for a gaussian in the dis-
%                        tance distribution, fraction of distance
%
%    max_fwhm          - maximum FWHM for a gaussian in the dis-
%                        tance distribution, fraction of distance
%                        range
%
%    max_mdep          - minimum DEER modulation depth
%
%    min_mdep          - maximum DEER modulation depth
%
%    expt              - background model, 'deer' or 'ridme'
%
%    max_brate         - maximum background signal decay
%                        rate, s^-1 (DEER backgrounds only)
%
%    min_brate         - minimum background signal decay 
%                        rate, s^-1 (DEER backgrounds only)
%
%    min_bdim          - minimum background dimensionality
%                        (DEER backgrounds only)
%
%    max_bdim          - maximum background dimensionality
%                        (DEER backgrounds only)
%
%    max_tshift        - maximum number of time discretisation
%                        points to shift the trace by, either
%                        forward or backward
%
% The function returns library.* structure with the following fields:
%
%    time_grid         - time grid (seconds) as a row vector
%
%    dist_grid         - distance grid (Angsrom) as a row vector
%
%    dist_distr_lib    - all distance distributions as a horizon-
%                        tal stack of column vectors
%
%    background_lib    - all background signals as a horizonal 
%                        stack of column vectors, shifted and 
%                        scaled to match DEER/RIDME traces
%
%    deer_noisy_lib    - all complete DEER/RIDME traces as a 
%                        horizonal stack of column vectors
%
%    deer_clean_lib    - DEER traces as they would come out, but
%                        without the noise; horizonal stack of 
%                        column vectors
% 
%    exchange_lib      - exchange interaction (MHz), a row vector
%                        containing the value for each example
%
%    parameters        - parameters array as received
%
% If a file name is provided, the structure is written into that file.
%
% Note: a strong NVidia GPU is required - at least a Titan V.
%
% jake.keeley@soton.ac.uk
% tajwar.choudhury@soton.ac.uk
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=deer_lib_gen.m>

function library=deer_lib_gen(file_name,parameters)
      
% Check consistency
grumble(file_name,parameters);

% Time grid for DEER traces (seconds)
time_grid=linspace(0,parameters.max_time,...
                     parameters.np_time);
dt=time_grid(2); parameters.time_grid=time_grid;

% Distance range supported by the current grid
[grid_rmin,grid_rmax]=dist_range(dt,parameters.max_time,parameters.expt);

% Distance grid for distributions (Angstroms)
dist_grid=linspace(grid_rmin,grid_rmax,parameters.np_dist);
parameters.dist_grid=dist_grid;

% Convert exchange couplings into rad/s
max_exch=parameters.max_exch*(pi/dt);
min_exch=parameters.min_exch*(pi/dt);
               
% Cached dipolar constants and Fresnel functions
cache_file_name=['cache_' md5_hash([time_grid dist_grid]) '.mat'];
if exist(cache_file_name,'file')
    load(cache_file_name,'D','FrC','FrS');
else
    D=zeros(size(dist_grid),'single');
    FrC=zeros(numel(time_grid),numel(dist_grid),'single');
    FrS=zeros(numel(time_grid),numel(dist_grid),'single');
    parfor n=1:numel(dist_grid)
        D(n)=xyz2dd([0 0 0],[0 0 dist_grid(n)],'E','E');
        FrC(:,n)=fresnelc(sqrt(6*D(n)*time_grid/pi))';
        FrS(:,n)=fresnels(sqrt(6*D(n)*time_grid/pi))';
    end
    save(cache_file_name,'D','FrC','FrS');
end

% Time grid to single precision
time_grid=single(time_grid)';

% Preallocate distance distributions
dist_distr_lib=zeros(parameters.np_dist,...
                     parameters.ntraces,'single');

% Preallocate background signals
background_lib=zeros(parameters.np_time,...
                     parameters.ntraces,'single');
                
% Preallocate DEER traces
deer_noisy_lib=zeros(parameters.np_time,...
                     parameters.ntraces,'single');
deer_clean_lib=zeros(parameters.np_time,...
                     parameters.ntraces,'single');
                 
% Preallocate the exchange coupling array
exchange_lib=zeros(1,parameters.ntraces,'single');

% Precompute DEER kernel components
if (min_exch==0)&&(max_exch==0)
    
    % Dipolar phase grid
    phase_grid=time_grid*D;
    
    % Zero exchange coupling: precompute full DEER kernel
    deer_kern=sqrt(pi./(6*phase_grid)).*(cos(phase_grid).*FrC+...
                                         sin(phase_grid).*FrS);
                                     
    % Clean-up and dummy variables for parfor
    clear('phase_grid','FrC','FrS','D'); 
    SQFrC=[]; SQFrS=[]; D=[];
                                      
else
    
    % Non-zero exchange: precompute what we can
    SQFrC=sqrt(pi./(6*time_grid*D)).*FrC;
    SQFrS=sqrt(pi./(6*time_grid*D)).*FrS;
    
    % Clean-up and dummy variables for parfor
    clear('FrC','FrS'); deer_kern=[];
    
end

% Save the kernel
parameters.kernel=deer_kern;
                 
% Generate distance distributions
parfor n=1:parameters.ntraces %#ok<*PFOUS,*PFBNS>
    
    % Random number of distance peaks
    ndists=randi(parameters.ndistmax);
    
    % Generate distance distribution
    dist_distr=zeros(parameters.np_dist,1);
    for k=1:ndists

        % Linear bias towards short distances
        trand=abs(rand('single')+rand('single')-1);

        % Randomly select a distance
        dist=grid_rmin+(grid_rmax-grid_rmin)*trand;
        
        % Get width range
        min_width=parameters.min_fwhm*dist;
        max_width=parameters.max_fwhm*(grid_rmax-grid_rmin);
        if min_width>max_width
            error('erroneous FWHM specification.');
        end
        
        % Linear bias towards narrow peaks
        trand=abs(rand('single')+rand('single')-1);
        
        % Randomly select a width
        fwhm=min_width+(max_width-min_width)*trand;
        sigma=fwhm/(2*sqrt(2*log(2)));
        
        % Add distance peak to distribution with a random amplitude
        dist_distr=dist_distr+rand('single')*normpdf(dist_grid,dist,sigma)';
        
    end
    
    % Normalise distance distribution
    dist_distr=dist_distr/sum(dist_distr);
    
    % Decide exchange coupling situation
    if (min_exch==0)&&(max_exch==0)
        
        % Use pre-computed kernel for zero exchange
        deer_ffact=deer_kern*dist_distr;
        
    else
    
        % Randomly select exchange coupling
        J=min_exch+(max_exch-min_exch)*rand(); exchange_lib(n)=J;
    
        % Use Kuprov's formula (http://dx.doi.org/10.1038/ncomms14842)
        deer_ffact=cos(time_grid*(D+J)).*SQFrC+...
                   sin(time_grid*(D+J)).*SQFrS;
        deer_ffact=deer_ffact*dist_distr;
        
    end
    
    % Fix the hole
    deer_ffact(1)=1;
                     
    % Randomly select modulation depth
    mdep=parameters.min_mdep+...
        (parameters.max_mdep-...
         parameters.min_mdep)*rand('single');
     
    if strcmp(parameters.expt,'deer')
             
        % Randomly select background dim
        bdim=parameters.min_bdim+...
            (parameters.max_bdim-...
             parameters.min_bdim)*rand();

        % Randomly select bg decay rate
        brate=parameters.min_brate+...
             (parameters.max_brate-...
              parameters.min_brate)*rand();

        % Build the background signal
        background=exp(-(brate*time_grid).^(bdim/3));
            
    elseif strcmp(parameters.expt,'ridme')
        
        % Randomly select parameters
        a=random('Normal',0,3);
        b=random('HalfNormal',max(-a/0.4,-a/2),3);

        % Get a scaled time axis
        scaled_time_grid=linspace(0,1,parameters.np_time)';

        % Build the background signal
        background=exp(-(a*scaled_time_grid+...
                         b*scaled_time_grid.^2));
                     
    else 
        
        % Complain and bomb out
        background=[]; %#ok<NASGU> 
        error('Unknown background model.');
            
    end
    
    % Apply the time grid shift
    if parameters.max_tshift>0

        % Randomly pick the shift amount
        shift_amount=randi(2*parameters.max_tshift+1)-parameters.max_tshift-1;

        % Apply forward shift
        if shift_amount>0
            deer_ffact=[flipud(deer_ffact(2:(shift_amount+1))); deer_ffact];
            deer_ffact=deer_ffact(1:(end-shift_amount));
            background=[flipud(background(2:(shift_amount+1))); background];
            background=background(1:(end-shift_amount));
        end

        % Apply backward shift
        if shift_amount<0
            shift_amount=abs(shift_amount);
            deer_ffact=[deer_ffact; single(lpredict(double(deer_ffact),...
                                           5*shift_amount,...
                                             shift_amount))];
            deer_ffact=deer_ffact((shift_amount+1):end);
            background=[background; single(lpredict(double(background),...
                                                    5*shift_amount,...
                                                      shift_amount))];
            background=background((shift_amount+1):end);
        end

    end

    % Mix background with form factor
    deer_clean=(1-mdep+mdep*deer_ffact).*background;
    
    % Match the modulation depth
    background=(1-mdep)*background;
        
    % Add a random amount of noise
    noise_line=rand('single')*parameters.noise_lvl*...
               mdep*randn(size(deer_clean),'single');
    deer_noisy=deer_clean+noise_line;
    
    % Box the noisy signal into [0 1]
    top_edge=max(deer_noisy); bot_edge=min(deer_noisy);
    deer_noisy=(deer_noisy-bot_edge)/(top_edge-bot_edge);
    deer_clean=(deer_clean-bot_edge)/(top_edge-bot_edge);
    background=(background-bot_edge)/(top_edge-bot_edge);
    
    % Normalize distance distribution
    dist_distr=dist_distr*parameters.np_dist;
    
    % Populate library arrays
    dist_distr_lib(:,n)=dist_distr;
    background_lib(:,n)=background;
    deer_noisy_lib(:,n)=deer_noisy;
    deer_clean_lib(:,n)=deer_clean;
    
end

% Build output data structure
library.parameters=parameters;
library.time_grid=time_grid;
library.dist_grid=dist_grid;
library.dist_distr_lib=dist_distr_lib;
library.deer_clean_lib=deer_clean_lib;
library.deer_noisy_lib=deer_noisy_lib;
library.background_lib=background_lib;
library.exchange_lib=exchange_lib;

% Save the data
if ~isempty(file_name)
    save(file_name,'library','-v7.3');
end

end

% Consistency enforcement
function grumble(file_name,parameters)
if ~isempty(file_name)&&~ischar(file_name)
    error('file_name must be a character string.');
end
if ~isfield(parameters,'max_time')
    error('DEER trace duration must be provided in parameters.max_time');
end
if ~isfield(parameters,'min_exch')
    error('minimum fractional exchange coupling must be provided in parameters.min_exch');
end
if ~isfield(parameters,'max_exch')
    error('maximum fractional exchange coupling must be provided in parameters.max_exch');
end
if ~isfield(parameters,'ntraces')
    error('number of traces must be provided in parameters.ntraces');
end
if ~isfield(parameters,'ndistmax')
    error('maximum number of distance peaks must be provided in parameters.ndistmax');
end
if ~isfield(parameters,'np_time')
    error('number of time digitisation points must be provided in parameters.np_time')
end
if ~isfield(parameters,'np_dist')
    error('number of distance digitisation points must be provided in parameters.np_dist');
end
if ~isfield(parameters,'noise_lvl')
    error('noise level must be provided in parameters.noise_lvl');
end
if ~isfield(parameters,'min_fwhm')
    error('minimum distance peak FWHM must be provided in parameters.min_fwhm');
end
if ~isfield(parameters,'max_fwhm')
    error('maximum distance peak FWHM must be provided in parameters.max_fwhm');
end
if ~isfield(parameters,'min_mdep')
    error('minimum modulation depth must be provided in parameters.min_mdep');
end
if ~isfield(parameters,'max_mdep')
    error('maximum modulation depth must be provided in parameters.max_mdep');
end
if ~isfield(parameters,'min_brate')
    error('minimum BG decay rate must be provided in parameters.min_brate');
end
if ~isfield(parameters,'max_brate')
    error('maximum BG decay rate must be provided in parameters.max_brate');
end
if ~isfield(parameters,'min_bdim')
    error('minimum BG dimensionality must be provided in parameters.min_bdim');
end
if ~isfield(parameters,'max_bdim')
    error('maximum BG dimensionality must be provided in parameters.max_bdim');
end
if ~isfield(parameters,'expt')
    error('experiment type must be provided in parameters.expt');
end
end

% "I know rejection is painful, and writing these emails has always 
%  been difficult for me - I also didn't get a similar prize as a PhD
%  student. A small and ironic solace may be in the fact that one day
%  you could be running the committee awarding them..."
%
% from IK's email template to unsuccessful 
% Bruker Thesis Prize applicants

