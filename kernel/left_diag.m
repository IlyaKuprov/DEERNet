% Generates a weight matrix descrambler for a particular layer in
% a neural network using maximum diagonality criterion described
% in (https://www.pnas.org/doi/10.1073/pnas.2016917118). Syntax:
%
%                P=left_diag(W,method,n_iter,guess)
%
% Parameters:
%
%    W       - layer weight matrix, must be square
%
%    method  - 'max_diag_sum' finds a transformation that creates
%              P*W with maximum diagonal sum; 'max_diag_normsq'
%              finds a transformation that creates P*W with maxi-
%              mum diagonal norm square.
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
%    P       - the matrix accomplishing the transformation when
%              it is multiplied into W from the left.
%   
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=left_diag.m>

function P=left_diag(W,method,n_iter,guess)

% Check consistency
grumble(W,method,n_iter);

% Decide problem dimensions
out_dim=size(W,1);
opt_dim=(out_dim^2-out_dim)/2;

% Lower triangle index array
lt_idx=tril(true(out_dim),-1);

% Default guess is zero
if ~exist('guess','var')
    guess=zeros(opt_dim,1);
else
    guess=guess(lt_idx);
end
   
% Precompute some steps
U=eye(out_dim);

% Regularisation signal
function [eta,eta_grad]=reg_sig(q)
    
    % Form the generator
    Q=zeros(out_dim); 
    Q(lt_idx)=q; Q=Q-Q';
    
    % Re-use the inverse
    iUpQ=inv(U+Q);
    
    % Run Cayley transform
    P=iUpQ*(U-Q); %#ok<MINV>
    
    % Precompute diag(P*W)
    DPW=diag(P*W);
    
    switch method
        
        case 'max_diag_sum'
    
            % Compute diagonal sum
            eta=-sum(DPW);
    
            % Compute diagonal sum gradient
            eta_grad=iUpQ'*(W')*(U+P)';
            
        case 'max_diag_normsq'
    
            % Compute diagonal norm-square
            eta=-norm(DPW,2)^2;
            
            % Compute diagonal norm-square gradient
            eta_grad=2*iUpQ'*(W'.*DPW)*(U+P)';
            
    end
    
    % Antisymmetrise the gradient
    eta_grad=eta_grad-transpose(eta_grad);
    
    % Extract the lower triangle
    eta_grad=eta_grad(lt_idx);

end

% Optimisation
options=optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
                     'MaxIterations',n_iter,'MaxFunctionEvaluations',inf,...
                     'FiniteDifferenceType','central',...
                     'SpecifyObjectiveGradient',true,'HessianApproximation',...
                     'finite-difference','SubproblemAlgorithm','cg');
q=fmincon(@reg_sig,guess,[],[],[],[],-inf(opt_dim,1),+inf(opt_dim,1),[],options); 

% Form descramble generator
Q=zeros(out_dim); Q(lt_idx)=q; Q=Q-Q';

% Run Cayley transform
P=(U+Q)\(U-Q);

end

% Consistency enforcement
function grumble(W,method,n_iter)
if (~isnumeric(W))||(~isreal(W))||...
   (size(W,1)~=size(W,2))
    error('W must be a square real matrix.');
end
if ~ischar(method)
    error('method must be a character string.');
end
if (~isnumeric(n_iter))||(~isreal(n_iter))||...
   (~isscalar(n_iter))||(mod(n_iter,1)~=0)||...
   (n_iter<1)
    error('n_iter must be a positive real integer.');
end
end

% "Animal courage was Lord Nelson's sole merit," said Lord Howe of his
% greatest lieutenant, "his private character was most disgraceful."
% Napoleon meanwhile had 27 mistresses, and David Lloyd George, when
% asked whether he was taking Mrs Lloyd George to the Paris Peace Con-
% ference, replied: "Would you take sandwiches to a banquet?"
%
% Andrew Roberts

