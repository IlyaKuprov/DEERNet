% Trains a network and writes a file with the specified num- 
% ber. If the specified file exists, the function will load
% that network file and continue training. Syntax:
%
%       train_one_net(parameters,batch_size,file_number)
%
% Parameters:
%
%    parameters  - library and network specification struc-
%                  ture (see the header of deer_lib_gen.m)
%
%    batch_size  - the number of training examples randomly
%                  generated at each iteration; must fit in
%                  memory
%
%    file_number - the network object will be saved into a
%                  file with this number as the name
%
% Outputs:
%
%    this function writes a .mat file with the network 
%    object and the parameters structure
%
% Note: a strong NVidia GPU is required - at least a Titan V
%
% Note: if the file pre-exists, the network will be loaded and
%       used as the initial guess.
%
% jake.keeley@soton.ac.uk
% tajwar.choudhury@soton.ac.uk
% s.g.worswick@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=train_one_net.m>

function train_one_net(parameters,batch_size,file_number)

% Check consistency
grumble(parameters,batch_size,file_number);

% Store training parameters
parameters.net_name=num2str(file_number);

% Set network file name
net_file_name=[num2str(file_number) '.mat'];

% Create a datastore instance
trainDS=DEERDatastore(parameters);
trainDS.MiniBatchSize=batch_size;

% Get the network
if isfile(net_file_name)
    
    % Try to load pre-existing net
    load(net_file_name,'net');
    if ~exist('net','var')
        error([net_file_name ' does not contain a network.']);
    else
        disp(['Initial guess loaded from ' net_file_name]);
    end
    layers=net.Layers;
    
else
    
    % Call network constructor
    if isfield(parameters,'npt_acq')&&...
       (parameters.npt_acq~=parameters.np_time)

        % Inform the user
        disp(['Training a (' num2str(parameters.npt_acq) ')-' ...
              num2str(parameters.np_time) '-'                 ...
              num2str(parameters.np_dist) ' DEERVet...']);
        
        % Get a blank DEERVet
        layers=dist_vet(parameters.np_dist,parameters.npt_acq);
        
    else

        % Inform the user
        disp(['Training a (' num2str(parameters.np_time) ')-' ...
              num2str(parameters.np_time) '-'                 ...
              num2str(parameters.np_dist) ' DEERNet...']);
        
        % Get a blank DEERNet
        layers=dist_net(parameters.np_time,parameters.np_dist);
        
    end
    
end

% Set training options
train_opts=trainingOptions('adam','MaxEpochs',50,'Verbose',true,'Shuffle','never',...
                           'InitialLearnrate',0.01,'LearnRateSchedule','piecewise',...
                           'LearnRateDropPeriod',5,'LearnRateDropFactor',0.5,...
                           'MiniBatchSize',batch_size,'ExecutionEnvironment','gpu',...
                           'L2Regularization',1e-6,'ResetInputNormalization',false,...
                           'Plots','none','DispatchInBackground',true);

% Train network
[net,train_stats]=trainNetwork(trainDS,layers,train_opts);

% Update parameter array
parameters.ntraces=0;
library=deer_lib_gen('',parameters);
parameters=library.parameters;
                        
% Save network
save([num2str(file_number) '.mat'],'net','train_stats','train_opts','parameters');

end

% Consistency enforcement
function grumble(parameters,batch_size,file_number)
if (~isnumeric(batch_size))||(~isreal(batch_size))||...
   (~isscalar(batch_size))||(batch_size<1000)||(mod(batch_size,1)~=0)
    error('batch_size must be a large positive integer.');
end
if (~isnumeric(file_number))||(~isreal(file_number))||...
   (~isscalar(file_number))||(file_number<=0)||...
   (mod(file_number,1)~=0)
    error('file_number must be a positive integer.');
end
if ~isfield(parameters,'training_input')
    error('training input type must be specified in parameters.training_input');
end
end

% Remember that rights are moral principles which define and 
% protect a man's freedom of action, but impose no obligati-
% ons on other men.
%
% Ayn Rand

