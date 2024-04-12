% Returns an untrained distance distribution DEERNet for pro-
% cessing sparsely sampled data. Syntax:
%
%               layers=dist_vet(npt_out,npt_acq)
%
% Parameters:
%
%     npt_out - the dimension of the output vector
%
%     npt_acq - the number of points acquired in the
%               sparsely sampled DEER trace; the input
%               vector dimension will be five times
%               this number to accommodate the lego
%
% Outputs:
%
%     layers  - an untrained network layout 
%
% i.kuprov@soton.ac.uk
% jake.keeley@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dist_vet.m>

function layers=dist_vet(npt_out,npt_acq)

% Check consistency
grumble(npt_out,npt_acq);

% Start with vector input layer
layers=featureInputLayer(5*npt_acq,'Normalization','none');
                   
% Six FC layers with softplus activations                   
layers=[layers; fullyConnectedLayer(5*npt_acq,'BiasLearnRateFactor',0,...
                                              'Bias',zeros(5*npt_acq,1));
                batchNormalizationLayer(); softplusLayer('Name','SP1')];
layers=[layers; fullyConnectedLayer(5*npt_acq,'BiasLearnRateFactor',0,...
                                              'Bias',zeros(5*npt_acq,1));
                batchNormalizationLayer(); softplusLayer('Name','SP2')];
layers=[layers; fullyConnectedLayer(5*npt_acq,'BiasLearnRateFactor',0,...
                                              'Bias',zeros(5*npt_acq,1));
                batchNormalizationLayer(); softplusLayer('Name','SP3')];            
layers=[layers; fullyConnectedLayer(5*npt_acq,'BiasLearnRateFactor',0,...
                                              'Bias',zeros(5*npt_acq,1));
                batchNormalizationLayer(); softplusLayer('Name','SP4')];            
layers=[layers; fullyConnectedLayer(5*npt_acq,'BiasLearnRateFactor',0,...
                                              'Bias',zeros(5*npt_acq,1));
                batchNormalizationLayer(); softplusLayer('Name','SP5')];            
layers=[layers; fullyConnectedLayer(npt_out,'BiasLearnRateFactor',0,...
                                            'Bias',zeros(npt_out,1));
                batchNormalizationLayer(); softplusLayer('Name','SP6')];            
            
% Probability renormalisation followed by regression
layers=[layers; renormLayer(); regressionLayer()];

end

% Consistency enforcement
function grumble(npt_out,npt_acq)
if (~isnumeric(npt_out))||(~isreal(npt_out))||...
   (~isscalar(npt_out))||(mod(npt_out,1)~=0)||(npt_out<1)
    error('npt_out must be a positive real integer.');
end
if (~isnumeric(npt_acq))||(~isreal(npt_acq))||...
   (~isscalar(npt_acq))||(mod(npt_acq,1)~=0)||(npt_acq<1)
    error('npt_acq must be a positive real integer.');
end
end

% Begin anywhere.
%
% John Cage, 
% a musician.

