% Converts Bruker Elexsys data sets into DEERNet input format. Adapted
% from DeerAnalysis package. The DEER signal is read from the .DTA file
% and the time axis is reconstructed using parameters read from the ac-
% companying .DSC file. The DEER signal is trimmed to ensure that the
% maximum occurs at t=0. Syntax:
% 
%     [deer_trace,time_axis]=elexsys2deernet(file,trunc,atrunc)
%
% Parameters:
% 
%         file  - filename for the .DSC and .DTA Elexsys 
%                 files, extension omitted
%
%         trunc - optional, a vector of two integers, specifying
%                 the number of points to be truncated from the 
%                 start and the end of the signal
%
%        atrunc - optional; if set to false, disables automatic
%                 left edge truncation at the maximum
%
% Outputs:
% 
%   deer_trace  - DEER trace, a column vector
%
%    time_axis  - time axis, a column vector (seconds)
%
% gunnar.jeschke@phys.chem.ethz.ch
% s.g.worswick@soton.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=elexsys2deernet.m>

function [deer_trace,time_axis]=elexsys2deernet(file,trunc,atrunc)

% Set default truncation
if ~exist('trunc','var'), trunc=[0 0]; end
if ~exist('atrunc','var'), atrunc=true(); end

% Check consistency
grumble(file,trunc,atrunc);

% Read the DTA file
fileID=fopen([file '.DTA'],'r','s'); disp(' ');
deer_trace=fread(fileID,inf,'float64'); fclose(fileID);
disp('========= DEERNet 2.9 data import from ELEXSYS =========');
disp(['DTA file: ' file '.DTA']);

% Read the time axis
if exist([file '.XGF'],'file')
    
    % Read XGF if present
    fileID=fopen([file '.XGF'],'r','s');
    time_axis=fread(fileID,inf,'float64'); 
    disp(['XGF file: ' file '.XGF']); fclose(fileID);
    
elseif exist([file '.DSC'],'file')

    % Otherwise read DSC
    fileID=fopen([file '.DSC'],'rt');
    DSC=textscan(fileID,'%s','Delimiter','\n');
    disp(['DSC file: ' file '.DSC']); fclose(fileID);

    % Extract the time axis from DSC file
    time_axis=linspace(get_value(DSC,'XMIN'),...
                       get_value(DSC,'XWID'),...
                       get_value(DSC,'XPTS'))';

else

    % Complain and bomb out
    error('.XGF or .DSC file with the same name as .DTA must be present.');

end

% Complex versus real data
if numel(deer_trace)==2*numel(time_axis)
    
    % Deinterleave real and imaginary parts
    deer_trace=reshape(deer_trace,[2 numel(time_axis)])';
    real_part=deer_trace(:,1); imag_part=deer_trace(:,2);
    
    % Apply truncation and time shift
    real_part=real_part((1+trunc(1)):(end-trunc(2)));
    imag_part=imag_part((1+trunc(1)):(end-trunc(2)));
    time_axis=time_axis((1+trunc(1)):(end-trunc(2)));
    time_axis=time_axis-time_axis(1);
    
    % Set optimiser options
    options=optimoptions('fminunc','FiniteDifferenceType','central',...
                         'Display','off','MaxIterations',1000,...
                         'MaxFunctionEvaluations',inf);
                     
    % Minimise squared norm of imaginary part
    scaling=norm(real_part+1i*imag_part,2);
    phi=fminunc(@(phi)targetfun(real_part/scaling,...
                                imag_part/scaling,phi),0,options);
    
    % Get new real and imaginary parts
    new_real=cos(phi)*real_part-sin(phi)*imag_part;
    new_imag=sin(phi)*real_part+cos(phi)*imag_part;
    
    % Inform the user
    disp(['Complex data autophased, phi = ' ...
          num2str(180*phi/pi) ' degrees:']);
    disp(['   |old_real| = ' num2str(norm(real_part,2)/scaling) ...
          ', |old_imag| = ' num2str(norm(imag_part,2)/scaling)]);
    disp(['   |new_real| = ' num2str(norm(new_real,2)/scaling) ...
          ', |new_imag| = ' num2str(norm(new_imag,2)/scaling)]);
      
    % Make sure real part is correct way up
    deer_trace=new_real*sign(sum(new_real));
        
else
    
    % When the data is real, only truncate and shift
    deer_trace=deer_trace((1+trunc(1)):(end-trunc(2)));
    time_axis=time_axis((1+trunc(1)):(end-trunc(2)));
    time_axis=time_axis-time_axis(1);

end

% Manual truncation diagnostics
if nargin>1

    % Report manual data truncation settings
    disp('Manual data truncation:');
    disp(['   ' num2str(trunc(1)) ' points deleted from the start,']);
    disp(['   ' num2str(trunc(2)) ' points deleted from the end,']);
    disp( '   time axis shifted to start at zero.');

end

% Automatic truncation
if atrunc

    % Locate signal maximum
    [~,index]=max(deer_trace);
    disp( 'Automatic data truncation:');
    disp(['   signal maximum found at point ' num2str(index) ',']);

    % Apply further truncation
    deer_trace=deer_trace(index:end);
    time_axis=time_axis(index:end);
    time_axis=time_axis-time_axis(1);

    % Inform the user
    if index>1, disp('   prior points deleted,'); end
    disp( '   time axis shifted to start at zero.');

end

% Convert time to seconds
time_axis=1e-9*time_axis;

% Close the section
disp('========================================================');

end

% Find parameter values in DSC file using Bruker mnemonics
function value=get_value(DSC,mnem)

% Find line with mnemonic
mask=~cellfun(@isempty,strfind(DSC{1},mnem));
line_str=DSC{1}(mask);

% Extract value from line string 
value=str2double(strrep(line_str{1},mnem,''));

end

% Autophasiong functional
function imnormsq=targetfun(real_part,imag_part,phi)

% Squared norm of the imaginary part
imnormsq=norm(sin(phi)*real_part+cos(phi)*imag_part,2)^2;

end

% Consistency enforcement
function grumble(file,trunc,atrunc)
if ~ischar(file), error('file must be a character string.'); end
if (~isnumeric(trunc))||(~isreal(trunc))||(numel(trunc)~=2)||...
   (any(trunc<0,'all'))||(any(mod(trunc,1)~=0,'all'))
    error('trunc must be a vector of two non-negative integers.');
end
if (~islogical(atrunc))||(~isscalar(atrunc))
    error('atrunc must be a logical scalar, use true() or false()');
end
end
            
% A good scientist is a person with original ideas. A good 
% engineer is a person who makes a design that works with 
% as few original ideas as possible. There are no prima 
% donnas in engineering.
% 
% Freeman Dyson

