% Gadolinium datasets from Akiva Feintuch and Daniella Goldfarb.
%
% Figure 5, https://doi.org/10.1039/c6cp00829a
%
% daniella.goldfarb@weizmann.ac.il
% akiva.feintuch@weizmann.ac.il

function example_set_gadolinium()

% Load the deta
load('./data_gadolinium/gd_deer.mat','gd_deer');

% Run the processing
deernet(gd_deer.deer_e,1e-6*gd_deer.t);
set(gcf,'Name','Offset: 1.09 GHz','NumberTitle','off');
drawnow();
deernet(gd_deer.deer_d,1e-6*gd_deer.t);
set(gcf,'Name','Offset: 0.747 GHz','NumberTitle','off');
drawnow();
deernet(gd_deer.deer_c,1e-6*gd_deer.t);
set(gcf,'Name','Offset: 0.469 GHz','NumberTitle','off');
drawnow();
deernet(gd_deer.deer_b,1e-6*gd_deer.t);
set(gcf,'Name','Offset: 0.363 GHz','NumberTitle','off');
drawnow();
deernet(gd_deer.deer_a,1e-6*gd_deer.t);
set(gcf,'Name','Offset: 0.106 GHz','NumberTitle','off');
drawnow();

end

