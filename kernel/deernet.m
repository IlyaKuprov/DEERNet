% Uses pre-trained neural networks to extract background signals, modula-
% tion depths, and distance distributions from DEER data. See our papers:
%
%                https://doi.org/10.1126/sciadv.aat5218
%                https://doi.org/10.1073/pnas.2016917118
%                https://doi.org/10.1016/j.jmr.2022.107186
%
% for further information. Syntax:
% 
%             dataset=deernet(input_traces,input_axis,options)
% 
% Required parameters:
% 
%    input_traces - experimental DEER trace(s), phased into pure absor-
%                   ption and cropped on the left to make the first po-
%                   int correspond to the echo modulation maximum; for
%                   sparsely sampled data, the missing points should be
%                   set to NaN. If multiple traces are supplied as col-
%                   umns of a matrix, they will be processed assuming
%                   common background dimension and decay rate.
%
%    input_axis   - experimental time axis in seconds, a column vector
%                   that must start at zero, have a uniform time step,
%                   and no NaN elements in it - even when the data is
%                   sparsely sampled
%
% Optional parameters (using new Matlab option syntax, literally as below): 
%
%    expt='ridme'       - requests RIDME type background, the default
%                         is DEER type background function
%
%    bg_dim_range=[a b] - constrains background signal dimension to 
%                         the interval [a b], default is [3.0,3.5],
%                         the training range was [2.0,3.5]; this op-
%                         tion only applies to expt='deer'.
%                      
%    do_jacobian=true   - requests Jacobian calculation; the run time
%                         will become much longer
%
% Output:
%
% either a figure (if there are no output parameters), or a data struc- 
% cture with the following fields:
%
%    input_traces   - DEER/RIDME trace(s), as supplied by the user
%
%    input_axis     - time axis, as supplied by the user
%
%    ntraces        - number of DEER/RIDME traces in the input
%
%    nsmpls         - number of non-NaN elements in the input DEER
%                     or RIDME trace
%
%    net_type       - 'net' for uniformly sampled input data, 'vet' 
%                     for sparsely sampled input data
%
%    resamp_axis    - time axis resampled to match DEERNet input dimen-
%                     sion, same as the input time axis for sparsely
%                     sampled inputs; this axis is used by background
%                     and retrocalculation outputs
%
%    resamp_traces  - DEER/RIDME trace(s) resampled to match DEERNet 
%                     input dimensions, same as the input trace for 
%                     sparsely sampled inputs
%
%    nnets          - number of networks in the netset that did the
%                     processing and the statistics
%
%    expt           - background model selected at input
%
%    bg_rates       - (DEER only) background decay rate returned by
%                     the retrofit of the data from each network
%
%    bg_dims        - (DEER only) background dimension returned by 
%                     the retrofit of the data from each network
%
%    backgs_av      - arithmetic mean of the background signal over the
%                     outputs of the networks in the current netset
%
%    backgs_lb      - 95% confidence interval lower bound for the back-
%                     ground signal obtained from netset statistics
%
%    backgs_ub      - 95% confidence interval upper bound for the back-
%                     ground signal obtained from netset statistics
%
%    retros_av      - arithmetic mean of the trace retrofit(s) over the
%                     outputs of the networks in the current netset
%    
%    retros_lb      - 95% confidence interval lower bound for the retrofit
%                     signal obtained from netset statistics
%
%    retros_ub      - 95% confidence interval upper bound for the retrofit
%                     signal obtained from netset statistics
%
%    mdpths_av      - arithmetic mean of the modulation depth(s) over the
%                     outputs of the networks in the current netset
%
%    mdpths_st      - standard deviation of modulation depth predictions
%                     obtained from netset statistics
%
%    dist_ax        - the distance axis matching the distances distribution
%                     and satisfying the constraints imposed by the timing
%                     of the input signal and the limits of DEERNet
%
%    dist_av        - arithmetic mean of the distance distributions(s) over
%                     the netset
%
%    dist_lb        - 95% confidence interval lower bound for the distance
%                     distribution(s) obtained from netset statistics
%
%    dist_ub        - 95% confidence interval upper bound for the distance
%                     distribution(s) obtained from netset statistics
%
% Note: for non-uniform sampling, the default netset can handle the 
%       following combinations of time grids and sampling schedules:
%
%          512-point time grid, 128 points sampled
%          512-point time grid, 64  points sampled
%
% jkeeley@mathworks.com
% tajwar.choudhury@soton.ac.uk
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=deernet.m>

function dataset=deernet(input_traces,input_axis,varargin)

% Enforce Matlab release
if isMATLABReleaseOlderThan('R2023a')
    error('DEERNet requires Matlab R2023a or later.');
end

% Enforce existential toolboxes
if ~exist([matlabroot '/toolbox/parallel'],'dir')
    error('DEERNet requires Parallel Computing Toolbox.');
end
if ~exist([matlabroot '/toolbox/nnet'],'dir')
    error('DEERNet requires Deep Learning Toolbox.');
end
if ~exist([matlabroot '/toolbox/rl'],'dir')
    error('DEERNet requires Reinforcement Learning Toolbox.');
end
if ~exist([matlabroot '/toolbox/optim'],'dir')
    error('DEERNet requires Optimisation Toolbox.');
end
if ~exist([matlabroot '/toolbox/stats'],'dir')
    error('DEERNet requires Statistics and Machine Learning Toolbox.');
end
if ~exist([matlabroot '/toolbox/map'],'dir')
    error('DEERNet requires Mapping Toolbox.');
end

% Banner and input parser object
disp(newline); p=inputParser; 
disp('===================== DEERNet 2.9 ======================');

% Default to DEER background model
addOptional(p,'expt','deer');

% Default to no Jacobian (expensive)
addOptional(p,'do_jacobian',false);

% Default to typical DEER BG dimensions
addOptional(p,'bg_dim_range',[3.0 3.5]);

% Parse the input or use the defaults
parse(p,varargin{:}); expt=p.Results.expt;
if strcmp(expt,'deer')
   
    % Report background type
    disp('Background model: DEER (stretched exponential)');

    % Report parameter ranges
    bg_dim_range=p.Results.bg_dim_range;
    if (bg_dim_range(1)==3.0)&&(bg_dim_range(2)==3.5)

        % The most common dimension range
        disp(['DEER BG dimension range: [' ...
               num2str(bg_dim_range(1)) ',' ...
               num2str(bg_dim_range(2)) '] (default)']);

    elseif isnumeric(bg_dim_range)&&isreal(bg_dim_range)&&...
           (numel(bg_dim_range)==2)&&(bg_dim_range(1)>=2.0)&&...
           (bg_dim_range(2)<=3.5) % Training database range

        % User-specified dimension range
        disp(['DEER BG dimension range: [' ...
               num2str(bg_dim_range(1)) ',' ...
               num2str(bg_dim_range(2)) '] (user-specified)']);

    else

        % Complain and bomb out
        error('incorrect DEER background dimension range specification.');

    end

elseif strcmp(expt,'ridme')

    % Report background type
    disp('Background model: RIDME (shifted gaussian)');

else

    % Complain and bomb out
    error('unknown experiment type.');

end

% Report the error analysis path
do_jacobian=p.Results.do_jacobian;
if do_jacobian, disp('Jacobian calculation enabled'); end

% Kernel set to dipolar, to be extended
if ~exist('netset','var'), netset='dd'; end
disp('Kernel type: dd (pure dipole-dipole)');

% Check consistency
grumble(input_traces,input_axis,expt,netset);

% Convert inputs to double precision and store
dataset.input_traces=double(input_traces);
dataset.input_axis=double(input_axis);

% Note the number of traces
ntraces=size(input_traces,2); dataset.ntraces=ntraces;
if ntraces>1
    disp(['Number of traces in the batch: ' ...
           num2str(ntraces) ' (shared background)']);
else
    disp(['Number of traces in the batch: ' ...
           num2str(ntraces)]);
end

% Check sampling schedules
nsmpls=unique(sum(~isnan(input_traces),1)); dataset.nsmpls=nsmpls;
if numel(nsmpls)~=1, error('sparse traces have different sample counts.'); end

% Select appropiate network type
if nsmpls~=numel(input_axis)
    
    % Report the sampling statistics
    net_type='vet'; disp(['Data sampling: sparse, npts=' ...
                           num2str(numel(input_axis))    ...
                          ', nsamp=' num2str(nsmpls)]);
    
    % First sample must be present
    if any(isnan(input_traces(1,:)),'all')
        error('with sparse sampling, the first element of the dataset cannot be missing.');
    end

    % Last sample must be present
    if any(isnan(input_traces(end,:)),'all')
        error('with sparse sampling, the last element of the dataset cannot be missing.');
    end

else

    % Report the sampling statistics
    net_type='net'; disp(['Data sampling: uniform, npts='       ...
                           num2str(numel(input_axis)) ', tmax=' ...
                           num2str(max(1e6*input_axis)) ' us']);

end
dataset.net_type=net_type;

% Locate the network ensemble parent directory
location=which(['library_' netset '.m']);
if isempty(location)
    error(['parameters not found, make sure library_' ...
            netset 'is on MATLAB path.'])
end
location=[location(1:(end-20)) filesep 'netset'];
location=[location filesep net_type '_distan_' netset];

% Specify the network directory naming pattern
patt=expt+"-("+digitsPattern+")-"+digitsPattern+"-512";

% Search the file system for network directories
dir_listing={dir(location).name}; 
match_list=dir_listing(matches(dir_listing,patt));
available_topols=reshape(str2double(extract(match_list,digitsPattern)),[],3);

% Choose the best network topology for the data
required_topol=[nsmpls numel(input_axis) NaN];
[val,idx]=min(sum(abs(available_topols-required_topol),2,'omitnan'));
if ~isempty(idx)
    disp(['Net layout, expt-(samp)-npts-nout: ' match_list{idx}]);
    best_topol=available_topols(idx,:); 
else
    error(['no suitable networks found in ' location]);
end

% Topology must be an exact match if data is sparse
if strcmp(net_type,'vet') && val~=0
    error(['exact topology match not found, ' ...
           'sparse data processing unavailable.'])
end

% Update the location
location=fullfile(location,match_list{idx});

% Generate the distance axes supported
% by the sampling conditions
switch net_type

    case 'net'

        % Resample input data to match network input dimension
        resamp_axis=linspace(input_axis(1),input_axis(end),best_topol(2))';
        resamp_traces=interp1(input_axis,input_traces,resamp_axis,'pchip');
        disp(['Data resampled from ' num2str(numel(input_axis)) ' to ' ...
               num2str(numel(resamp_axis)) ' time points']);

        % Box each signal into [0 1] to match training
        top_edges=max(resamp_traces,[],1); bot_edges=min(resamp_traces,[],1);
        net_inputs=(resamp_traces-bot_edges)./(top_edges-bot_edges);

        % Get the distance range supported by the network dimensions
        [deernet_rmin,deernet_rmax]=dist_range(resamp_axis(2),...
                                               resamp_axis(end),expt);

        % Generate the distance axis supported by the networks
        dist_axis=linspace(deernet_rmin,deernet_rmax,best_topol(end))';

        % Get the distance range supported by the input sampling
        [input_rmin,input_rmax]=dist_range(dataset.input_axis(2),...
                                           dataset.input_axis(end),expt);

        % Generate distance range truncation mask
        dist_keep_mask=(dist_axis>=input_rmin)&(dist_axis<=input_rmax);

    case 'vet'

        % Resampling is impossible with missing points
        resamp_axis=dataset.input_axis;
        resamp_traces=dataset.input_traces;

        % Box each signal into [0 1] to match training
        top_edges=max(resamp_traces,[],1,'omitnan'); 
        bot_edges=min(resamp_traces,[],1,'omitnan');
        boxed_traces=(resamp_traces-bot_edges)./(top_edges-bot_edges);

        % Get the distance range supported by the uniform time grid
        [uniform_rmin,uniform_rmax]=dist_range(resamp_axis(2),...
                                               resamp_axis(end),expt);

        % Generate the distance axis supported by the networks
        dist_axis=linspace(uniform_rmin,uniform_rmax,best_topol(end))';

        % Preallocate loop outputs
        net_inputs=zeros(nsmpls.*5,ntraces);
        sparse_rmin=zeros(1,ntraces); sparse_rmax=zeros(1,ntraces);

        % Get the distance range supported by the input sampling
        for n=1:ntraces

            % Extract the schedule
            schedule=find(~isnan(boxed_traces(:,n)));

            % Extract the data points
            data_pts=boxed_traces(schedule,n);

            % Prepare the network input
            mask_lin=2*(schedule/best_topol(2)-0.5);
            mask_sqr=4*(schedule/best_topol(2)-0.5).^2-0.5;
            net_inputs(:,n)=[mask_lin; mask_sqr; data_pts;
                             data_pts.*mask_lin; data_pts.*mask_sqr];

            % Get the distance range supported by the smapling
            sample_times=resamp_axis(schedule);
            sparse_dt=min(diff(sample_times)); sparse_tmax=max(sample_times);
            [sparse_rmin(n),sparse_rmax(n)]=dist_range(sparse_dt,sparse_tmax,expt);

        end

        % Generate distance range truncation mask
        dist_keep_mask=(dist_axis>max(sparse_rmin))&...
                       (dist_axis<min(sparse_rmax));

end

% Report distance extents
disp(['Visible distance range: [' ...
       num2str(min(dist_axis(dist_keep_mask))) ',' ...
       num2str(max(dist_axis(dist_keep_mask))) '] Angstrom']);

% Keep resampled axis and data
dataset.resamp_axis=resamp_axis;
dataset.resamp_traces=resamp_traces;

% Get scaling factors for preconditioning
scaling_factors=max(resamp_traces,[],1,'omitnan');

% Set upper and lower bounds for background parameters
switch expt

    case 'deer'
        
        % DEER background has physical motivation
        lb=[-inf(ntraces,1); -inf(ntraces,1);   0; bg_dim_range(1)];
        ub=[+inf(ntraces,1); +inf(ntraces,1); inf; bg_dim_range(2)];

    case 'ridme'

        % RIDME background is more libertarian
        lb=[-inf(ntraces,1); -inf(ntraces,1); -inf; -inf];
        ub=[+inf(ntraces,1); +inf(ntraces,1); +inf; +inf];

end

% Count networks in the ensemble
net_files=dir([location filesep '*.mat']); 
nnets=numel(net_files); dataset.nnets=nnets;

% Create finite difference increments and preallocate Jacobian stack
if do_jacobian
    dm=sqrt(eps('single'))*eye(size(net_inputs,1),'single'); dx=2*sqrt(eps('single'));
    jacob=zeros(best_topol(end),size(net_inputs,1),ntraces,nnets,'single');
else
    dm=[]; dx=[]; jacob=zeros([0 0 ntraces nnets]); % Appease parfor syntax checker
end

% Preallocate network output sets
backgs=zeros(best_topol(2),ntraces,nnets);
distds=zeros(best_topol(3),ntraces,nnets);
retros=zeros(best_topol(2),ntraces,nnets);
mdpths=zeros(1,ntraces,nnets);
bg_param_a=zeros(nnets,1);
bg_param_b=zeros(nnets,1);

% Loop over the networks
parfor n=1:nnets

    % Get the distance distribution network file name
    dist_net_name=[net_files(n).folder filesep net_files(n).name];

    % Call the distance distribution network
    distd=process_using(net_inputs,dist_net_name);

    % Compute Jacobians
    if do_jacobian
        for k=1:ntraces
            jacob(:,:,k,n)=(process_using(net_inputs(:,k)+dm,dist_net_name)-...
                            process_using(net_inputs(:,k)-dm,dist_net_name))/dx;
            jacob(:,:,k,n)=jacob(:,:,k,n)/scaling_factors(k); %#ok<PFBNS>
        end
    end

    % Catch softplus over- and underflows
    if any(isnan(distd),'all')
        disp(['Location: ' dist_net_name]);
        error('this network returns NaN, check your input data.')
    end

    % Load the kernel from the network file
    kernel=getfield(getfield(matfile(dist_net_name),'parameters'),'kernel');

    % Compute the form factors and remove singularities
    ffs=kernel*(distd./sum(distd,1)); ffs(1,:)=1; ffs=double(ffs);

    % Scale input data for retrofit
    scaled_traces=resamp_traces./scaling_factors;

    % Set the optimisation target functional
    lsq_err=@(x)double(sum(nonnans(scaled_traces-signal_model(x,ffs,ntraces,expt)).^2,'all'));

    % Set the initial guess
    guess=[0.50*(1+rand())*ones(ntraces,1);
           0.25*(1+rand())*ones(ntraces,1);
           0.50*(1+rand()); 2.0+1.5*rand()];

    % Set optimiser options
    options=optimoptions('fmincon','Algorithm','interior-point',...
                         'FiniteDifferenceType','central',...
                         'HessianApproximation','lbfgs',...
                         'InitTrustRegionRadius',0.01,...
                         'ScaleProblem',false,'Display','off',...
                         'MaxIterations',1000,...
                         'MaxFunctionEvaluations',inf);

    % Run the optimisation
    x=fmincon(lsq_err,guess,[],[],[],[],lb,ub,[],options);

    % Caclulate the retrofits and backgrounds
    retros(:,:,n)=scaling_factors.*signal_model(x,ffs,ntraces,expt);
    backgs(:,:,n)=scaling_factors.*signal_model(x,0*ffs,ntraces,expt);

    % Extract the modulation depths
    mdpths(1,:,n)=x((ntraces+1):(2*ntraces));

    % Weigh the distance distributions
    distds(:,:,n)=mdpths(1,:,n).*distd;

    % Keep background parameters
    bg_param_a(n)=x(end-1); bg_param_b(n)=x(end);

end

% Store experiment type
dataset.expt=expt;

% Store BG parameters
switch expt

    case 'deer'
        
        % Store background decay 
        % rates and dimensions
        dataset.bg_rates=bg_param_a; 
        dataset.bg_dims=bg_param_b;

    case 'ridme'

        % Store RIDME background 
        % a and b parameters
        dataset.bg_ridme_a=bg_param_a; 
        dataset.bg_ridme_b=bg_param_b;

end

% Confidence intervals on backgrounds at 95%
dataset.backgs=backgs; dataset.backgs_av=mean(backgs,3);
dataset.backgs_lb=mean(backgs,3)-2*std(backgs,[],3);
dataset.backgs_ub=mean(backgs,3)+2*std(backgs,[],3);

% Confidence intervals on retrocalcs at 95%
dataset.retros=retros; dataset.retros_av=mean(retros,3);
dataset.retros_lb=mean(retros,3)-2*std(retros,[],3);
dataset.retros_ub=mean(retros,3)+2*std(retros,[],3);

% Standard deviation calculation for modulation depths
dataset.mdpths=mdpths; dataset.mdpths_av=mean(mdpths,3);
dataset.mdpths_st=std(mdpths,[],3);

% Confidence intervals on distance distribution at 95%
dataset.distds=distds; dataset.dist_av=mean(distds,3);
dataset.dist_lb=mean(distds,3)-2*std(distds,[],3);
dataset.dist_ub=mean(distds,3)+2*std(distds,[],3);

% Cut the bottom off the lower bound
dataset.dist_lb(dataset.dist_lb<0)=0;

% Apply distance range truncation
dataset.dist_ax=dist_axis(dist_keep_mask);
dataset.distds=dataset.distds(dist_keep_mask,:,:);
dataset.dist_av=dataset.dist_av(dist_keep_mask,:);
dataset.dist_lb=dataset.dist_lb(dist_keep_mask,:);
dataset.dist_ub=dataset.dist_ub(dist_keep_mask,:);
if do_jacobian
    dataset.jacobian=jacob(dist_keep_mask,:,:,:);
end

% Run quality control
dataset=quality_control(dataset);

% Do the plotting only if no output
if nargout==0, deerplot(dataset); end

end

% Consistency enforcement
function grumble(input_traces,input_axis,expt,netset)
if (~ischar(expt))||(~ismember(expt,{'deer','ridme'}))
    error('expt must be ''deer'' or ''ridme''.');
end
if (~ischar(netset))||(~ismember(netset,{'dd'}))
    error('only netset=''dd'' is supported at the moment.');
end
if (~isnumeric(input_traces))||(~isreal(input_traces))||...
   (~ismatrix(input_traces))
    error('input traces must be a real column vector or matrix');
end
if (~isnumeric(input_axis))||(~isreal(input_axis))||...
   (~iscolumn(input_axis))
    error('time axis must be a real column vector.');
end
if (numel(input_axis)~=size(input_traces,1))
    error('dimensions of input_axis and input_traces must agree')
end
if norm(diff(input_axis,2),1)/norm(diff(input_axis,1),1)>1e-3
    error('time axis must have a uniform step.');
end
if input_axis(1)~=0, error('time axis must start at zero.'); end
if max(input_axis)>1e-4, error('time axis units must be seconds.'); end
end

% To succeed as a band you only need one good song, 
% and to succeed as a scientist - one good paper. 
% Look at Aerosmith and Al Redfield.
%
% Anonymous colleague at 
% the ENC Conference

