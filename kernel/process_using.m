% Applies the specified trained neural network file to the supplied 
% DEER trace. Syntax:
%
%          answers=process_using(deer_traces,net_file_name)
%
% Parameters:
%
%   deer_traces  -  DEER trace(s) as a column vector or a 
%                   matrix with multiple columns. The num-
%                   ber of rows must match the input size 
%                   of the neural network.
%
%   net_file_name - the name of the file containing a neu-
%                   ral network object, including full
%                   path and the .mat extension.
%
% Outputs:
%
%   answers       - neural network output(s) as a column 
%                   vector or a matrix with multiple col-
%                   umns. The number of rows matches the 
%                   output size of the neural network.
%
% jkeeley@mathworks.com
% i.kuprov@soton.ac.uk
% s.g.worswick@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=process_using.m>

function answers=process_using(deer_traces,net_file_name)

% Check consistency
grumble(deer_traces,net_file_name);

% Try to load the network
load(net_file_name,'net');
if ~exist('net','var')
    error([net_file_name ' does not contain a neural network.']);
end

% Run through the specified net on CPU
answers=predict(net,deer_traces','MiniBatchSize',size(deer_traces,2),...
                                 'ExecutionEnvironment','cpu');

% Reshape answers to match the input
answers=transpose(squeeze(answers));

end

% Consistency enforcement
function grumble(deer_traces,net_file_name)
if (~isnumeric(deer_traces))||...
   (~isreal(deer_traces))||(~ismatrix(deer_traces))
    error('deer_traces must be a real column vector or matrix.');
end
net_file_name=convertStringsToChars(net_file_name);
if (~ischar(net_file_name))||...
   (~strcmp(net_file_name((end-3):end),'.mat'))
    error('net_file_name must be a character string ending in .mat');
end
end

% Those who say that all cultures are equal never explain why 
% the results of those cultures are so grossly unequal.
%
% Thomas Sowell

