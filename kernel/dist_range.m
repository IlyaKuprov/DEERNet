% Returns the mathematical bounds for the distance range covered
% by a DEER trace with the specified number of points and durati-
% on, assuming uniform time step. Syntax:
%
%              [rmin,rmax]=dist_range(dt,tmax,expt)
%
% Parameters:
%
%      dt  -  interval duration of the uniform time grid
%
%    tmax  -  the time (seconds) of the last point in the
%             grid if the first point is zero
%
%    expt  -  'ridme' or 'deer'
%
% Outputs:
%
%    rmin  -  recommended minimum distance for the distance
%             grid, Angstrom 
%
%    rmax  -  recommended maximum distance for the distance
%             grid, Angstrom 
%
% i.kuprov@soton.ac.uk
% jkeeley@mathworks.com
%
% <https://spindynamics.org/wiki/index.php?title=dist_range.m>

function [rmin,rmax]=dist_range(dt,tmax,expt)

% Check consistency
grumble(dt,tmax);

% Fundamental constants
hbar=1.054571628e-34;
mu0=4*pi*1e-7; 
gamma=-1.76085963023e11;

% Decide the sampling
switch expt
    
    case 'deer' % Nyquist-Shannon minimal sampling
        
        % Only half dipolar period needs to be present
        % because Pake pattern has twice dipolar frequency
        rmax=1e10*(2*(mu0/(4*pi))*(gamma^2*hbar/(2*pi))*tmax)^(1/3);

        % A minimum of four intervals per dipolar period
        % because Pake pattern has twice dipolar frequency
        rmin=1e10*(4*(mu0/(4*pi))*(gamma^2*hbar/(2*pi))*dt)^(1/3);
    
    case 'ridme'   % Empirically sensible sampling
        
        % Require the presence of a full dipolar period
        rmax=1e10*((mu0/(4*pi))*(gamma^2*hbar/(2*pi))*tmax)^(1/3);

        % A minimum of eight intervals per dipolar period
        rmin=1e10*(8*(mu0/(4*pi))*(gamma^2*hbar/(2*pi))*dt)^(1/3);

    otherwise
        
        % Complain and bomb out
        error('unknown safety level.');

end
end

% Consistency enforcement
function grumble(dt,tmax)
if (~isnumeric(dt))||(~isreal(dt))||...
   (~isscalar(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
if (~isnumeric(tmax))||(~isreal(tmax))||...
   (~isscalar(tmax))||(tmax<=0)
    error('tmax must be a positive real scalar.');
end
end

% Now I challenge people to tell me a joke that's not of- 
% fensive and I can find something offensive in it. 
%   - "Why did the chicken cross the road?" 
%   - "How dare you, my chicken died yesterday!"
%
% Ricky Gervais

