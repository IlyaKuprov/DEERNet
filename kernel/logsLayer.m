% Logsigmoidal activation function layer. Syntax:
%
%                     layer=logsLayer()
%
% Outputs:
%
%    layer - Matlab Deep Learning Toolbox layer object
%
% jkeeley@mathworks.com
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=logsLayer.m>

classdef logsLayer < nnet.layer.Layer
    
    methods
        
        % Constructor
        function layer=logsLayer()
            layer.Description="Logsigmoidal function";
            layer.Type="Logsig";
        end
        
        % Predictor
        function Z=predict(layer,X) %#ok<*INUSL>
            Z=1./(1+exp(-X));
        end

    end
    
end

% "In any case, to wear an improper expression on your face (to
%  look incredulous when a victory was announced, for example)
%  was itself a punishable offence. There was even a word for it
%  in Newspeak: facecrime, it was called." 
%
% George Orwell, 1984

