mpopt = mpoption('opf.init_from_mpc', 1, 'opf.flow_lim', 'P');
runmyopf('case5_custom', mpopt)