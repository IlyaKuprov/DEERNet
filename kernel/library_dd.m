% Training database parameters for the purely dipolar DEERNet. Max 
% time and max distance are in a sliding relationship, so the trace
% length is notional - picked to enable intuition about parameters.
% Exchange couplings are specified as fractions of the maximum freq-
% uency that is representable on the time grid. Syntax:
%
%                      parameters=library_dd()
%
% Outputs:
%
%    parameters - data structure required by deer_lib_gen.m 
%                 and other DEERNet design functions
%
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=library_dd.m>

function parameters=library_dd()

% Empirical parameters
parameters.max_time=2e-6;    % DEER trace length (notional)
parameters.ndistmax=4;       % Maximum number of distance peaks
parameters.min_exch=0.0;     % Minimum fractional exchange coupling
parameters.max_exch=0.0;     % Maximum fractional exchange coupling
parameters.min_bdim=2.0;     % Minimum background dimensionality
parameters.max_bdim=3.5;     % Maximum background dimensionality
parameters.min_fwhm=0.025;   % Minimum FWHM, fraction of distance
parameters.max_fwhm=0.500;   % Maximum FWHM, fraction of distance range
parameters.min_mdep=0.01;    % Minimum modulation depth
parameters.max_mdep=1.00;    % Maximum modulation depth
parameters.noise_lvl=0.05;   % Maximum RMS noise, fraction of MD
parameters.min_brate=0.0e6;  % Minimum background decay rate, s^-1
parameters.max_brate=0.5e6;  % Maximum background decay rate, s^-1
parameters.max_tshift=0;     % Maximum time shift, points

% Diagnostic messages
disp('Network set: purely dipolar modulation, arbitrary distance distribution.');

end

% Soviet cosmonaut Sergei Krikalev was in space when the Soviet 
% Union fell apart in 1991. Unable to return home, he ended up
% having to stay in space until further notice. He eventually
% came back to Earth on 25 March 1992, after 10 months in orbit,
% to a nation that was very different to what it had been when
% he had left. The Soviet Union had fractured into 15 countries,
% and even his hometown of Leningrad had become St. Petersburg.
% Interestingly, at the time, Krikalev was supposed to serve in
% the military reserves, and was almost issued a warrant for de-
% sertion - before the army realised that their reserve soldier
% was not even on the planet.
%
% Ran Levi, "The Awful and Wonderful History 
%            of the Mir Space Station"

