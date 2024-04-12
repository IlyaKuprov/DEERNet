% Returns the non-NaN elements of a numerical array. Syntax:
%
%                       nums=nonnans(a)
%
% Parameters:
%
%     a    - a numeric array of any dimension
%
% Output:
%
%     nums - non-NaN elements of the input array,
%            a column vector
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=nonnans.m>

function nums=nonnans(a)

% Check consistency
grumble(a);

% Simple indexing
nums=a(~isnan(a));

end

% Consistency enforcement
function grumble(a)
if ~isnumeric(a)
    error('the input must be numeric.');
end
end

% Ты ебалась со стайером,
% А потом со спринтером.
% А я ебусь с таймером
% И ещё с принтером.
% 
% Женя Лесин

