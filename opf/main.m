clear;
mpopt = mpoption('opf.init_from_mpc', 1, 'opf.flow_lim', 'P');
runmyopf('case5_custom', mpopt)
% mpopt = mpoption('opf.init_from_mpc', 1, 'opf.flow_lim', 'P');
% runmyopf('case39', mpopt)