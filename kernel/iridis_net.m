% Trains one DEERNet on an Iridis5 supercomputer node at the
% University of Southampton. Syntax:
%
%                 iridis_net(n,batch_size)
%
% Parameters:
%
%     n  - network file number. If the file pre-exists, it
%          will be loaded and used as an initial guess
%
%     batch_size - batch size; 2^13 for rough training,
%                  2^17 for fine training
%
% Outputs:
%
%     the function writes a file containing a trained net
%
% Notes: this function will not run on any other system 
%        without appropriate modifications; use it as a
%        starting point.
%
% Notes: reduce parameters.npt_acq to train a DEERVet; set
%        parameters.expt='ridme' to train RIDME networks.
%
% i.kuprov@soton.ac.uk
%
% #NGRUM #NWIKI

function iridis_net(n,batch_size)

% Make sure jobs don't collide
mdcs_dir=['/scratch/ik1r11/net_p' num2str(n)];
mkdir(mdcs_dir); c=parcluster('local');
c.JobStorageLocation=mdcs_dir; parpool(c,20);

% Get dataset parameters
parameters=library_dd();

% Set background model
parameters.expt='deer';

% Number of elements in the data vector
parameters.np_time=512;

% Number of elements acquired
parameters.npt_acq=512;

% Number of elements in the prediction
parameters.np_dist=512;

% Run training against noisy data
parameters.training_input='noisy_data';
train_one_net(parameters,batch_size,n);

end

% As for the review experience, remember that story about 
% Pavlov's dogs? Conditional reflexes and things. There's
% a less well known experiment when dogs are punished and
% rewarded in a way that's uncorrelated with what they do.
% The dogs eventually develop schizophrenia... that's how
% that project review in Brussels made me feel.
%
% (from IK's email to the project team,
%  after the final review meeting on an
%  EU grant, in which the outcomes were
%  praised by the European Commission)

