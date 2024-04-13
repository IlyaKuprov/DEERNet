% DEER signal model used for DEERNet and DEERVet backcalcu-
% lations. Includes the form factor, the modulation depth,
% and the background. Syntax:
%
%          sig=signal_model(x,ffs,ntraces,expt)
%
% Parameters:
%
%    x       - array of parameters: N multiplicati-
%              ve prefactros, N modulation depths,
%              one background decay rate, one back-
%              ground dimension
%
%    ffs     - DEER form factors for each trace in
%              the common-background batch, a mat-
%              rix with individual form factors in
%              columns
%
%    ntraces - number of DEER traces in the common-
%              background batch
%
%    expt    - 'deer' for DEER background, 'ridme'
%              for RIDME background model
%
% Outputs:
%
%    sig     - DEER traces of the common-background
%              batch, a matrix with individual DEER
%              traces in columns
%
% jkeeley@mathworks.com
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=signal_model.m>

function sig=signal_model(x,ffs,ntraces,expt)

% Check consistency
grumble(x,ffs,ntraces,expt);

switch expt

    case 'deer'

        % Individual prefactors
        alphas=x(1:ntraces)';

        % Individual modulation depths
        mus=x((ntraces+1):(2*ntraces))';

        % Common decay rate
        k=x(2*ntraces+1);

        % Common BG dimension
        bgd=x(2*ntraces+2);

        % Compute signal array
        scaled_time=linspace(0,1,size(ffs,1))';
        sig=alphas.*((1-mus)+mus.*ffs).*...
        exp(-(k*scaled_time).^(bgd/3));
    
    case 'ridme'
        
        % Individual prefactors
        alphas=x(1:ntraces)';

        % Individual modulation depths
        mus=x((ntraces+1):(2*ntraces))';

        % Common a parameter
        a=x(2*ntraces+1);

        % Common b parameter 
        b=x(2*ntraces+2);

        % Compute signal array
        scaled_time=linspace(0,1,size(ffs,1))';
        sig=alphas.*((1-mus)+mus.*ffs).*...
        exp(-(a*scaled_time+b*scaled_time.^2));
    
end

end

% Consistency enforcement
function grumble(x,ffs,ntraces,expt)
if (~ismember(expt,{'deer','ridme'}))
    error('expt must be ''deer'' or ''ridme''.')
end
if (~isnumeric(ntraces))||(~isreal(ntraces))||...
   (~isscalar(ntraces))||(ntraces<1)||(mod(ntraces,1)~=0)
    error('ntraces must be a positive integer scalar.');
end
if (~isnumeric(x))||(~isreal(x))
    error('x must be a real numeric array.');
end
if (~isnumeric(ffs))||(~isreal(ffs))
    error('ffs must be a real numeric array.');
end
end

% Charlie Chaplin was one day at a fair in the United States, where
% a principal attraction was a competition as to who could best imi-
% tate the Charlie Chaplin walk. The real Charlie Chaplin thought 
% there might be a chance for him so he entered for the performance
% [...] He was a frightful failure and came in twentieth.
%
% Mary Pickford

