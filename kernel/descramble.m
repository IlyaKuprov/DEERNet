% Generates a weight matrix descrambler for a particular layer in
% a neural network using Tikhonov smoothness criterion. The parti-
% culars are described in 
%
%        https://www.pnas.org/doi/10.1073/pnas.2016917118
%
% Syntax:
%
%                     P=descramble(S,n_iter,guess)
%
% Parameters:
%
%    S       - a matrix containing, in its columns, the outputs
%              of the preceding layers of the neural network for
%              a (preferably large) number of reasonable inputs
%
%    n_iter  - maximum number of Newton-Raphson interations, 400
%              is generally sufficient
%
%    guess   - [optional] the initial guess for the descrambling
%              transform generator (lower triangle is used), a
%              reasonable choice is a zero matrix (default)
%
% Output:
%
%    P       - descrambling matrix. In the case when the network 
%              is wiretapped before the activation function, i.e.
%
%                             S = Wf(W...f(Wf(WX)))
%
%              matrix P descrambles the output dimension of the
%              left-most W. In the case when the network is wire-
%              tapped after the activation function, i.e.
%
%                            S = f(Wf(W...f(Wf(WX))))
%
%              matrix inv(P) descrambles the input dimension of
%              the weight matrix of the subsequent layer.
%   
% i.kuprov@soton.ac.uk
% j.amey@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=descramble.m>

function [P,Q]=descramble(S,n_iter,guess)

% Check consistency
grumble(S,n_iter);

% Decide problem dimensions
out_dim=size(S,1);
opt_dim=(out_dim^2-out_dim)/2;

% Lower triangle index array
lt_idx=tril(true(out_dim),-1);

% Default guess is zero
if ~exist('guess','var')
    guess=zeros(opt_dim,1);
else
    guess=guess(lt_idx);
end
   
% Get second derivative operator
[~,D]=fourdif(out_dim,2);

% Precompute some steps, scale, and move to GPU
SST=gpuArray(S*S'); SST=out_dim*SST/norm(SST,2);
DTD=gpuArray(D'*D); DTD=out_dim*DTD/norm(DTD,2);
U=eye([out_dim out_dim],'gpuArray');

% Regularisation signal
function [eta,eta_grad]=reg_sig(q)
    
    % Form the generator
    Q=zeros([out_dim out_dim],'gpuArray'); 
    Q(lt_idx)=q; Q=Q-Q';
    
    % Re-use the inverse
    iUpQ=inv(U+Q);
    
    % Run Cayley transform
    P=iUpQ*(U-Q); %#ok<MINV>
    
    % Re-use triple product
    DTDPSST=DTD*P*SST;
    
    % Compute Tikhonov norm
    eta=trace(DTDPSST*P');
    
    % Compute Tikhonov norm gradient
    eta_grad=-2*iUpQ'*DTDPSST*(U+P)';
    
    % Antisymmetrise the gradient
    eta_grad=eta_grad-transpose(eta_grad);
    
    % Extract the lower triangle
    eta_grad=eta_grad(lt_idx);
    
    % Move back to CPU
    eta=gather(eta); 
    eta_grad=gather(eta_grad);

end

% Optimisation
options=optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
                     'MaxIterations',n_iter,'MaxFunctionEvaluations',inf,...
                     'FiniteDifferenceType','central',...
                     'SpecifyObjectiveGradient',true,'HessianApproximation',{'lbfgs',5});
q=fmincon(@reg_sig,guess,[],[],[],[],-inf(opt_dim,1),+inf(opt_dim,1),[],options);

% Form descramble generator
Q=zeros(out_dim); Q(lt_idx)=q; Q=Q-Q';

% Run Cayley transform
P=(U+Q)\(U-Q);

end

% Consistency enforcement
function grumble(S,n_iter)
if (~isnumeric(S))||(~isreal(S))
    error('S must be a real matrix.');
end
if (~isnumeric(n_iter))||(~isreal(n_iter))||...
   (~isscalar(n_iter))||(mod(n_iter,1)~=0)||...
   (n_iter<1)
    error('n_iter must be a positive real integer.');
end
end

% =======================
% Dear ###,
% 
% you may have been wondering, for some of your past publications, 
% why the peer review process was taking such a long time. I would
% like to point out - if I may - that now you know.
%
% Best wishes,
% Ilya.
% =======================
%
% IK's reminder email to junior
% scientists who are dragging
% their feet on a paper review 

