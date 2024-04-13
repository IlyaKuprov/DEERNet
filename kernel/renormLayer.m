% Renormalisation layer for probability distributions. Makes 
% the mean value across the input vector equal to 1. The mean
% is used instead of the sum to stay around the sweet spot of
% single-precision arithmetic. Syntax:
%
%                    layer=renormLayer()
%
% Outputs:
%
%    layer - Matlab Deep Learning Toolbox layer object
%
% i.kuprov@soton.ac.uk
% jkeeley@mathworks.com
%
% <https://spindynamics.org/wiki/index.php?title=renormLayer.m>

classdef renormLayer < nnet.layer.Layer
    
    methods
        
        % Constructor
        function layer=renormLayer()
            layer.Description='Prob. Renorm. Layer';
            layer.Type='Renormalisation';
            layer.Name='Renorm';
            layer.NumInputs=1;
        end
        
        % Predictor
        function X=predict(layer,X) %#ok<*INUSL>
            X=X./mean(X,1);
        end

    end
    
end

% To a mathematician, a beautiful proof is one so quickly stated
% that the afternoon remains open for playing with Rubik's cubes
% and making fun of social scientists.
%
% SMBC Comics

