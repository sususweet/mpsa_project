function [results, success, raw] = myopf_execute(om, mpopt)
%OPF_EXECUTE  Executes the OPF specified by an OPF model object.
%   [RESULTS, SUCCESS, RAW] = OPF_EXECUTE(OM, MPOPT)
%
%   RESULTS are returned with internal indexing, all equipment
%   in-service, etc.
%
%   See also OPF, OPF_SETUP.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  setup  -----
%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
sdp = strcmp(alg, 'SDPOPF');

%% build user-defined costs
om = build_cost_params(om);

%% get indexing
[vv, ll, nn] = get_idx(om);

if mpopt.verbose > 0
    v = mpver('all');
    fprintf('\nMYOPF Version %s, %s', v.Version, v.Date);
end

%%-----  run DC OPF solver  -----
if dc
  if mpopt.verbose > 0
    fprintf(' -- DC Optimal Power Flow\n');
  end
  [results, success, raw] = dcopf_solver(om, mpopt);
else
  %%-----  run AC OPF solver  -----
  if mpopt.verbose > 0
    fprintf(' -- AC Optimal Power Flow\n');
  end

  %% if opf.ac.solver not set, choose best available option
  if strcmp(alg, 'DEFAULT')
    if have_fcn('pdipmopf')
      alg = 'PDIPM';            %% PDIPM
    else
      alg = 'MYOPF';             %% MYOPF
    end
    mpopt = mpoption(mpopt, 'opf.ac.solver', alg);
  end

  %% ZIP loads?
  if (~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
          any(mpopt.exp.sys_wide_zip_loads.pw(2:3))) || ...
          (~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
          any(mpopt.exp.sys_wide_zip_loads.qw(2:3)))
    switch alg
    case {'PDIPM', 'TRALM', 'MINOPF', 'SDPOPF'}
      warning('opf_execute: ''%s'' solver does not support ZIP load model. Converting to constant power loads.', alg)
      mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', ...
                        struct('pw', [], 'qw', []));
    end
  end

  %% run specific AC OPF solver
  switch alg
    case 'MYOPF'
      [results, success, raw] = myopf_solver(om, mpopt);
    case 'MIPS'
      [results, success, raw] = mipsopf_solver(om, mpopt);
    case 'IPOPT'
      if ~have_fcn('ipopt')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see http://www.coin-or.org/projects/Ipopt.xml)', alg);
      end
      [results, success, raw] = ipoptopf_solver(om, mpopt);
    case 'PDIPM'
      if mpopt.pdipm.step_control
        if ~have_fcn('scpdipmopf')
          error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires SCPDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
        end
      else
        if ~have_fcn('pdipmopf')
          error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires PDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
        end
      end
      [results, success, raw] = tspopf_solver(om, mpopt);
    case 'TRALM'
      if ~have_fcn('tralmopf')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires TRALM (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
      [results, success, raw] = tspopf_solver(om, mpopt);
    case 'MINOPF'
      if ~have_fcn('minopf')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires MINOPF (see http://www.pserc.cornell.edu/minopf/)', alg);
      end
      [results, success, raw] = mopf_solver(om, mpopt);
    case 'FMINCON'
      if ~have_fcn('fmincon')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires FMINCON (Optimization Toolbox 2.x or later)', alg);
      end
      [results, success, raw] = fmincopf_solver(om, mpopt);
    case 'KNITRO'
      if ~have_fcn('knitro')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires KNITRO (see http://www.ziena.com/)', alg);
      end
      [results, success, raw] = ktropf_solver(om, mpopt);
    case 'SDPOPF'
      if ~have_fcn('yalmip')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires YALMIP (see http://users.isy.liu.se/johanl/yalmip/)', alg);
      end
      [results, success, raw] = sdpopf_solver(om, mpopt);
    otherwise
      error('opf_execute: MPOPT.opf.ac.solver = ''%s'' is not a valid AC OPF solver selection', alg);
  end
end
if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
    raw.output.alg = alg;
end

if success
  if ~dc && ~sdp
    %% copy bus voltages back to gen matrix
    results.gen(:, VG) = results.bus(results.gen(:, GEN_BUS), VM);

    %% gen PQ capability curve multipliers
    if ll.N.PQh > 0 || ll.N.PQl > 0
      mu_PQh = results.mu.lin.l(ll.i1.PQh:ll.iN.PQh) - results.mu.lin.u(ll.i1.PQh:ll.iN.PQh);
      mu_PQl = results.mu.lin.l(ll.i1.PQl:ll.iN.PQl) - results.mu.lin.u(ll.i1.PQl:ll.iN.PQl);
      Apqdata = userdata(om, 'Apqdata');
      results.gen = update_mupq(results.baseMVA, results.gen, mu_PQh, mu_PQl, Apqdata);
    end

    %% compute g, dg, f, df, d2f if requested by opf.return_raw_der = 1
    if mpopt.opf.return_raw_der
      %% move from results to raw if using v4.0 of MINOPF or TSPOPF
      if isfield(results, 'dg')
        raw.dg = results.dg;
        raw.g = results.g;
      end
      %% compute g, dg, unless already done by post-v4.0 MINOPF or TSPOPF
      if ~isfield(raw, 'dg')
        mpc = get_mpc(om);
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
        [g, geq, dg, dgeq] = opf_consfcn(results.x, om, Ybus, Yf, Yt, mpopt);
        raw.g = [ geq; g];
        raw.dg = [ dgeq'; dg'];   %% true Jacobian organization
      end
      %% compute df, d2f
      [f, df, d2f] = opf_costfcn(results.x, om);
      raw.df = df;
      raw.d2f = d2f;
    end
  end

  %% delete g and dg fieldsfrom results if using v4.0 of MINOPF or TSPOPF
  if isfield(results, 'dg')
    rmfield(results, 'dg');
    rmfield(results, 'g');
  end

  %% angle limit constraint multipliers
  if ~sdp && ll.N.ang > 0
    iang = userdata(om, 'iang');
    results.branch(iang, MU_ANGMIN) = results.mu.lin.l(ll.i1.ang:ll.iN.ang) * pi/180;
    results.branch(iang, MU_ANGMAX) = results.mu.lin.u(ll.i1.ang:ll.iN.ang) * pi/180;
  end
else
  %% assign empty g, dg, f, df, d2f if requested by opf.return_raw_der = 1
  if ~dc && mpopt.opf.return_raw_der
    raw.dg = [];
    raw.g = [];
    raw.df = [];
    raw.d2f = [];
  end
end

if ~sdp
  %% assign values and limit shadow prices for variables
  om_var_order = get(om, 'var', 'order');
  for k = 1:length(om_var_order)
    name = om_var_order(k).name;
    if getN(om, 'var', name)
      idx = vv.i1.(name):vv.iN.(name);
      results.var.val.(name) = results.x(idx);
    end
  end

  %% assign shadow prices for linear constraints
  om_lin_order = get(om, 'lin', 'order');
  for k = 1:length(om_lin_order)
    name = om_lin_order(k).name;
    if getN(om, 'lin', name)
      idx = ll.i1.(name):ll.iN.(name);
      results.lin.mu.l.(name) = results.mu.lin.l(idx);
      results.lin.mu.u.(name) = results.mu.lin.u(idx);
    end
  end


  %% assign values for components of user cost
  om_cost_order = get(om, 'cost', 'order');
  for k = 1:length(om_cost_order)
    name = om_cost_order(k).name;
    if getN(om, 'cost', name)
      results.cost.(name) = compute_cost(om, results.x, name);
    end
  end

  %% if single-block PWL costs were converted to POLY, insert dummy y into x
  %% Note: The "y" portion of x will be nonsense, but everything should at
  %%       least be in the expected locations.
  pwl1 = userdata(om, 'pwl1');
  if ~isempty(pwl1) && ~strcmp(alg, 'TRALM') && ~(strcmp(alg, 'PDIPM') && mpopt.pdipm.step_control)
    %% get indexing
    vv = get_idx(om);
    if dc
      nx = vv.iN.Pg;
    else
      nx = vv.iN.Qg;
    end
    y = zeros(length(pwl1), 1);
    raw.xr = [ raw.xr(1:nx); y; raw.xr(nx+1:end)];
    results.x = [ results.x(1:nx); y; results.x(nx+1:end)];
  end
end