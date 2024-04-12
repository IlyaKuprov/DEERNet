% Elimination of orientation selection by pulse position
% averaging as described in 
%
%      http://dx.doi.org/10.1016/j.jmr.2012.11.028
%
% Data kindly provided by Daniella Goldfarb.
%
% i.kuprov@soton.ac.uk

function example_set_orisel()

load('data_orisel/orisel_data.mat','deer_data');
deernet(mean([deer_data.(2) deer_data.(3) deer_data.(4)],2),deer_data.(1));
subplot(1,2,1); plot(1e6*deer_data.(1),deer_data.(2),'Color',[0.75 0.5 0.5]);
                plot(1e6*deer_data.(1),deer_data.(3),'Color',[0.75 0.5 0.5]);
                plot(1e6*deer_data.(1),deer_data.(4),'Color',[0.75 0.5 0.5]);
legend('hide'); drawnow();

end

