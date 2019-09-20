function mpc = case5_custom
mpc.version='2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA=100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1   1   160 80  0   0   1   1   0     100   1   1.1 0.9;
    2   1   200 100 0   0   1   1   0     100   1   1.1 0.9;
    3   1   370 130 0   0   1   1   0     100   1   1.1 0.9;
    4   2   0   0   0   0   1   1.05   0    1   1   1.1 0.9;
    5   3   0   0   0   0   1   1.05   0    1   1   1.1 0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
    4   500 0   300   -300   1.05    100 1   800 100   0   0   0   0   0   0   0   0   0   0   0;
    5   0   0   500   -210   1.05    100 1   800 100   0   0   0   0   0   0   0   0   0   0   0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    2   1   0.04    0.25    0.5     200   0   0   0       0   1   -360    360;
    3   1   0.1     0.35    0.0     65    0   0   0       0   1   -360    360;
    2   3   0.08    0.3     0.5     200   0   0   0       0   1   -360    360;
    3   5   0.0     0.03    0.0     500   0   0   1.05    0   1   -360    360;
    2   4   0.0     0.015   0.0     600   0   0   1.05    0   1   -360    360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
    2   0   0   3   0.00504395   2.004335   1200.6485;
    2   0   0   3   0.020055     5.00746    1857.201;
];