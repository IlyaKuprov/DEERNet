% Returns an untrained distance distribution DEERNet for pro-
% cessing fully sampled data. Syntax:
%
%                 layers=dist_net(np_in,np_out)
%
% Parameters:
%
%     np_in   - dimension of the input vector
%
%     np_out  - dimension of the output vector
%
% Outputs:
%
%     layers  - an untrained network layout 
%
% i.kuprov@soton.ac.uk
% jake.keeley@soton.ac.uk
% tajwar.choudhury@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dist_net.m>

function layers=dist_net(np_in,np_out)

% Check consistency
grumble(np_in,np_out);

% Start with vector input layer
layers=featureInputLayer(np_in,'Normalization','none');

% Six FC layers with batch norm and softplus activation
layers=[layers; fullyConnectedLayer(np_out,'Bias',zeros(np_out,1),...
                                           'BiasLearnRateFactor',0);
                batchNormalizationLayer(); softplusLayer('Name','SP1')];
layers=[layers; fullyConnectedLayer(np_out,'Bias',zeros(np_out,1),...
                                           'BiasLearnRateFactor',0);
                batchNormalizationLayer(); softplusLayer('Name','SP2')];
layers=[layers; fullyConnectedLayer(np_out,'Bias',zeros(np_out,1),...
                                           'BiasLearnRateFactor',0);
                batchNormalizationLayer(); softplusLayer('Name','SP3')];
layers=[layers; fullyConnectedLayer(np_out,'Bias',zeros(np_out,1),...
                                           'BiasLearnRateFactor',0);
                batchNormalizationLayer(); softplusLayer('Name','SP4')];
layers=[layers; fullyConnectedLayer(np_out,'Bias',zeros(np_out,1),...
                                           'BiasLearnRateFactor',0);
                batchNormalizationLayer(); softplusLayer('Name','SP5')];
layers=[layers; fullyConnectedLayer(np_out,'Bias',zeros(np_out,1),...
                                           'BiasLearnRateFactor',0);
                batchNormalizationLayer(); softplusLayer('Name','SP6')];

% Probability renormalisation followed by regression
layers=[layers; renormLayer(); regressionLayer()];

end

% Consistency enforcement
function grumble(np_in,np_out)
if (~isnumeric(np_in))||(~isreal(np_in))||...
   (~isscalar(np_in))||(mod(np_in,1)~=0)||(np_in<1)
    error('np_in must be a positive real integer.');
end
if (~isnumeric(np_out))||(~isreal(np_out))||...
   (~isscalar(np_out))||(mod(np_out,1)~=0)||(np_out<1)
    error('np_out must be a positive real integer.');
end
end

% What one programmer can do in one month, two
% programmers can do in two months.
%
% Fred Brooks

