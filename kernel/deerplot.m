% Plotting subsystem of DEERNet. Syntax:
%
%                  deerplot(dataset)
%
% Parameters:
%
%    dataset - DEERNet output data structure
%
% Outputs:
%
%    this function creates a figure; the best way to
%    save it for publication is to use exportgraphics() 
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=deerplot.m>

function deerplot(dataset)

% Check consistency
grumble(dataset);

% Count the traces
ntraces=size(dataset.input_traces,2);

% Primary data plotting style
if any(isnan(dataset.input_traces),'all')
    line_style='.';
else
    line_style='-';
end

% Loop over the traces
for n=1:ntraces
    
    % Make decisions
    plot_background=isfield(dataset,'resamp_axis')&&(~any(isnan(dataset.resamp_axis)))&&...
                    isfield(dataset,'backgs_av')&&(~any(isnan(dataset.backgs_av(:,n))))&&...
                    isfield(dataset,'backgs_lb')&&(~any(isnan(dataset.backgs_lb(:,n))))&&...
                    isfield(dataset,'backgs_ub')&&(~any(isnan(dataset.backgs_ub(:,n))));
    plot_mdepth=isfield(dataset,'mdpths_av')&&(~isnan(dataset.mdpths_av(n)))&&...
                isfield(dataset,'mdpths_st')&&(~isnan(dataset.mdpths_st(n)));

    % Get the figure going
    figure(); loc=get(gcf,'Position');
    old_fig_cent=[loc(1)+loc(3)/2 loc(2)+loc(4)/2];
    loc=get(0,'defaultfigureposition');
    old_fig_size=[loc(3) loc(4)]; new_fig_cent=old_fig_cent;
    new_fig_size=[1.5 0.75].*old_fig_size;
    set(gcf,'Position',[new_fig_cent-new_fig_size/2 new_fig_size]);
    set(gcf,'defaultTextInterpreter','LaTeX')
    set(gcf,'defaultLegendInterpreter','LaTeX')
    set(gcf,'defaultAxesTickLabelInterpreter','LaTeX')

    % Start the left panel
    subplot(1,2,1); hold on; box on;

    % Start the legends
    legend_left={}; legend_right={};

    % Plot the input data
    if isfield(dataset,'input_axis')&&...
       isfield(dataset,'input_traces')
        plot(1e6*dataset.input_axis,...
                 dataset.input_traces(:,n),...
                 line_style,'Color',[0.5 0.5 0.5]);
        legend_left{end+1}='echo modulation'; %#ok<AGROW>
    end

    % Plot the background
    if plot_background       
        plot(1e6*dataset.resamp_axis,...
                 dataset.backgs_av(:,n),'Color',[0 0 0.75]);
        legend_left{end+1}='background, mean'; %#ok<AGROW>
        patch([1e6*dataset.resamp_axis;    flip(1e6*dataset.resamp_axis)],...
              [    dataset.backgs_lb(:,n); flip(    dataset.backgs_ub(:,n))],...
              [0 0 0.75],'FaceAlpha',0.15,'EdgeColor','none');
        legend_left{end+1}='background, 95\%';  %#ok<AGROW>
    end
    
    % Display background parameters
    if plot_mdepth

        % Full background analysis
        disp(['Trace ' num2str(n) ', modulation depth, average: ' ...
                       num2str(dataset.mdpths_av(n))]);
        disp(['Trace ' num2str(n) ', modulation depth, st.dev.: ' ...
                       num2str(dataset.mdpths_st(n))]);
        disp(['Trace ' num2str(n) ', (init. ampl.)*MD, average: ' ...
                       num2str(dataset.input_traces(1,n)*dataset.mdpths_av(1,n))]);
        disp(['Trace ' num2str(n) ', (init. ampl.)*MD, st.dev.: ' ...
                       num2str(dataset.input_traces(1,n)*dataset.mdpths_st(1,n))]);
        if isfield(dataset,'bg_rates')
            disp(['Trace ' num2str(n) ', background decay rate (MHz), average: ' ...
                           num2str(mean(dataset.bg_rates))]);
            disp(['Trace ' num2str(n) ', background decay rate (MHz), st.dev.: ' ...
                           num2str(std(dataset.bg_rates))]);
        end
        if isfield(dataset,'bg_dims')
            disp(['Trace ' num2str(n) ', background dimension, average:  ' ...
                           num2str(mean(dataset.bg_dims))]);
            disp(['Trace ' num2str(n) ', background dimension, st.dev.:  ' ...
                           num2str(std(dataset.bg_dims))]);
        end
        if isfield(dataset,'bg_ridme_a')
            disp(['Trace ' num2str(n) ', RIDME BG parameter A, average:  ' ...
                           num2str(mean(dataset.bg_ridme_a))]);
            disp(['Trace ' num2str(n) ', RIDME BG parameter A, st.dev.:  ' ...
                           num2str(std(dataset.bg_ridme_a))]);
        end
        if isfield(dataset,'bg_ridme_b')
            disp(['Trace ' num2str(n) ', RIDME BG parameter B, average:  ' ...
                           num2str(mean(dataset.bg_ridme_b))]);
            disp(['Trace ' num2str(n) ', RIDME BG parameter B, st.dev.:  ' ...
                           num2str(std(dataset.bg_ridme_b))]);
        end
        disp(' '); % Blank line separator

    else

        % Refuse to output background data 
        disp(['Trace ' num2str(n) ' background quantification failed.']);

    end

    % Plot the retrofitted signal
    if isfield(dataset,'resamp_axis')&&...
       isfield(dataset,'retros_av')
        plot(1e6*dataset.resamp_axis,...
                 dataset.retros_av(:,n),'Color',[0 0.75 0]);
        legend_left{end+1}='back calc., mean'; %#ok<AGROW>
    end

    % Plot the retrofit confidence band
    if isfield(dataset,'resamp_axis')&&...
       isfield(dataset,'retros_lb')&&...
       isfield(dataset,'retros_ub')
         patch([1e6*dataset.resamp_axis;    flip(1e6*dataset.resamp_axis)],...
               [    dataset.retros_lb(:,n); flip(    dataset.retros_ub(:,n))],...
               [0 0.75 0],'FaceAlpha',0.15,'EdgeColor','none');
         legend_left{end+1}='back calc., 95\%'; %#ok<AGROW>
    end

    % Residual cosmetics on the left panel
    axis tight; set(gca,'Layer','top'); 
    xlabel('time, $\mu$s'); ylabel('amplitude, a.u.');
    legend(legend_left); grid on; 
    set(gca,'GridAlpha',1,'Layer','bottom',...
            'GridColor',[0.85 0.85 0.85]);

    % Start the right panel
    subplot(1,2,2); hold on; box on; title(''); 

    % Plot the distance distribution
    if isfield(dataset,'dist_ax')&&...
       isfield(dataset,'dist_av')&&...
       ~isnan(dataset.dist_av(1,n))
        plot(dataset.dist_ax,dataset.dist_av(:,n),'Color',[0 0 0.75]);
        legend_right{end+1}='mean'; %#ok<AGROW>
    end

    % Plot distance confidence bands
    if isfield(dataset,'dist_ax')&&...
       isfield(dataset,'dist_lb')&&...
       isfield(dataset,'dist_ub')&&...
       ~isnan(dataset.dist_lb(1,n))&&...
       ~isnan(dataset.dist_ub(1,n))
        patch([dataset.dist_ax;      flip(dataset.dist_ax)],...
              [dataset.dist_lb(:,n); flip(dataset.dist_ub(:,n))],...
              [0 0 0.75],'FaceAlpha',0.15,'EdgeColor','none');
       legend_right{end+1}='95\%'; %#ok<AGROW>
    end

    % Residual cosmetics on the right panel
    axis tight; set(gca,'Layer','top');
    xlabel('distance, $\rm{\AA}$');
    ylabel('$\mu{p(r)}$'); grid on; 
    set(gca,'GridAlpha',1,'Layer','bottom',...
            'GridColor',[0.85 0.85 0.85]);
    if ~isempty(legend_right)
        legend(legend_right);
    else
        xlim([0 50]); ylim([0 1]);
        text(25,0.5,'you must be joking',...
             'HorizontalAlignment','center');
    end

end

end

% Consistency enforcement
function grumble(dataset)
if ~isstruct(dataset)
    error('dataset must be a structure from deernet() or deervet()');
end
if (~isfield(dataset,'input_traces'))||(~isfield(dataset,'input_axis'))
    error('dataset.input_axis and dataset.input_traces must be present');
end
end

% 27 And when Jesus departed thence, two blind men followed Him, 
%    crying, and saying, Thou son of David, have mercy on us.
% 
% 28 And when He was come into the house, the blind men came to 
%    Him: and Jesus saith unto them, Believe ye that I am able
%    to do this? They said unto Him, Yea, Lord.
% 
% 29 Then touched He their eyes, saying, According to your faith
%    be it unto you.
% 
% Matthew 9

