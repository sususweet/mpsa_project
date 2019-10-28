function [results, success, raw] = myopf_solver(om, mpopt)
%SUSUSWEETOPF_SOLVER  Solves AC optimal power flow using Interior Point method provided in Modern Power System Analysis.
%       min F(X)
%        X
%
%   subject to
%
%       G(X) = 0            (nonlinear equalities)
%       H(X) <= 0           (nonlinear inequalities)
%       L <= A*X <= U       (linear constraints)
%       XMIN <= X <= XMAX   (variable bounds)
%
%   [RESULTS, SUCCESS, RAW] = MIPSOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options struct.
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   RESULTS is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .nln
%               .l  lower bounds on nonlinear constraints
%               .u  upper bounds on nonlinear constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .output solver specific output information
%
%   See also OPF, MIPS.

%   MATPOWER
%   Copyright (c) 2000-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% options
opt = mpopt.mips;
opt.verbose = mpopt.verbose;
if opt.feastol == 0
    opt.feastol = mpopt.opf.violation;  %% = MPOPT.opf.violation by default
end
if ~isfield(opt, 'cost_mult') || isempty(opt.cost_mult)
    opt.cost_mult = 1e-4;
end

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs
ng = size(gen, 1);
%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

%% linear constraints
[A, l, u] = linear_constraints(om);

%% bounds on optimization vars
% 系统中所有变量的约束条件，排列顺序为theta1~n,V1~n,P1~g,Q1~g
[x0, xmin, xmax] = getv(om);
% 修改x0的排列顺序，使之符合课本的变量顺序
% x0_P_tmp = x0(2*nb+1:2*nb+ng);
% x0_Q_tmp = x0(2*nb+1+ng:end);
x0_P_tmp = (gen(:, PMAX) + gen(:, PMIN)) ./ 2 ./ baseMVA;
x0_Q_tmp = (gen(:, QMAX) + gen(:, QMIN)) ./ 2 ./ baseMVA;
% x0_P_tmp = (gen(:, PMAX)) ./ baseMVA;
% x0_Q_tmp = (gen(:, QMAX)) ./ baseMVA;
x0(2*nb+1:end) = [];
% 系统中所有变量的约束条件，排列顺序为P1~g,Q1~g,theta1~n,V1~n
x0 = [x0_P_tmp; x0_Q_tmp; x0];

nvariable = size(x0, 1);

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point
% if mpopt.opf.init_from_mpc ~= 1
%     ll = xmin; uu = xmax;
%     ll(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
%     uu(xmax ==  Inf) =  1e10;
%     x0 = (ll + uu) / 2;         %% set x0 mid-way between bounds
%     k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
%     x0(k) = xmax(k) - 1;                    %% set just below upper bound
%     k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
%     x0(k) = xmin(k) + 1;                    %% set just above lower bound
%     Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);
%     x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
%     if ny > 0
%         ipwl = find(gencost(:, MODEL) == PW_LINEAR);
%         %     PQ = [gen(:, PMAX); gen(:, QMAX)];
%         %     c = totcost(gencost(ipwl, :), PQ(ipwl));
%         c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
%         x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
%         %     x0(vv.i1.y:vv.iN.y) = c + 0.1 * abs(c);
%     end
% end

%% ====================================================================
%% -----  run myopf  -----
%% Calculate Gap
% build initial variables
r_num = nvariable - nb + nl2;
r_num2 = r_num + nl2;
lslack = ones(r_num2, 1);
uslack = ones(r_num2, 1);
zLagr = ones(r_num2, 1);
wLagr = -0.5 * ones(r_num2, 1);
yLagr = 1e-10 * [ones(nb, 1); -ones(nb, 1)];
% lslack = 0.8*ones(r_num2, 1);
% uslack = 1.1*ones(r_num2, 1);
% zLagr = ones(r_num2, 1);
% wLagr = -1.5 * ones(r_num2, 1);
% yLagr = 1e-3 * [ones(nb, 1); -ones(nb, 1)];

iterCount = 1;
Gap = lslack' * zLagr - uslack' * wLagr;
Gap_arr(iterCount) = Gap;

%% ====================================================================
% Cg: nb, ng generator connection matrix
% (i; j)th element is 1 if generator j is located at bus i, 0 otherwise
Cg = zeros(nb, ng);
gen_point = gen(:, GEN_BUS);
for i = 1:size(gen_point,1)
    Cg(gen_point(i, 1), i) = -1;
end

while Gap >= 1e-6
    Gap_arr(iterCount) = Gap;
    miu = 0.1 * Gap / (2 * (r_num2));
    
    %% get variables from x0
    x0_Pg = x0(1:ng);
    x0_Qg = x0(ng + 1:2 * ng);
    x0_theta = x0(2 * ng + 1: 2 * ng + nb);
    x0_V = x0(2 * ng + nb + 1 : end);
    x0_V_theta = x0_V .* exp(1j * x0_theta);
    Ibus = Ybus * x0_V_theta;
    diagIbus = diag(Ibus);
    
    %% update result variables
    gen(:, PG) = x0_Pg * baseMVA;
    gen(:, QG) = x0_Qg * baseMVA;
    gen(:, VG) = x0(2*ng + nb + gen(:, GEN_BUS));
    bus(:, VA) = rad2deg(x0_theta);
    bus(:, VM) = x0_V;

    %% compute branch flows
    Sf = x0_V_theta(branch(:, F_BUS)) .* conj(Yf * x0_V_theta);  %% cplx pwr at "from" bus, p.u.
    St = x0_V_theta(branch(:, T_BUS)) .* conj(Yt * x0_V_theta);  %% cplx pwr at "to" bus, p.u.
    branch(:, PF) = real(Sf) * baseMVA;
    branch(:, QF) = imag(Sf) * baseMVA;
    branch(:, PT) = real(St) * baseMVA;
    branch(:, QT) = imag(St) * baseMVA;
    
    %% ====================================================================
    %% form function h first derivatives
    x0_V_theta_diag = diag(x0_V_theta);
    x0_V_norm_diag = diag(x0_V_theta ./ abs(x0_V_theta));
    dSbus_dVm = x0_V_theta_diag * conj(Ybus * x0_V_norm_diag) + conj(diagIbus) * x0_V_norm_diag;
    dSbus_dVa = 1j * x0_V_theta_diag * conj(diagIbus - Ybus * x0_V_theta_diag);
    
    h_grad_matrix = -[Cg', zeros(ng, nb);
        zeros(ng, nb), Cg';
        real(dSbus_dVa.'), imag(dSbus_dVa.');
        real(dSbus_dVm.'), imag(dSbus_dVm.')];
    

    % % define
    % f = branch(:, F_BUS);       %% list of "from" buses
    % t = branch(:, T_BUS);       %% list of "to" buses
    %
    % % compute currents
    % If = Yf * x0_V_theta;
    % It = Yt * x0_V_theta;
    %
    % diagVf      = diag(x0_V_theta(f));
    % diagIf      = diag(If);
    % diagVt      = diag(x0_V_theta(t));
    % diagIt      = diag(It);
    % diagV       = diag(x0_V_theta);
    % diagVnorm   = diag(x0_V_norm);
    % temp1       = zeros(nl, nb);
    % temp1(sub2ind([nl,nb], (1:nl)', f)) = x0_V_theta(f);
    % temp2       = zeros(nl, nb);
    % temp2(sub2ind([nl,nb], (1:nl)', f)) = x0_V_norm(f);
    % temp3       = zeros(nl, nb);
    % temp3(sub2ind([nl,nb], (1:nl)', t)) = x0_V_theta(t);
    % temp4       = zeros(nl, nb);
    % temp4(sub2ind([nl,nb], (1:nl)', t)) = x0_V_norm(t);
    % dSf_dVa = 1j * (conj(diagIf) * temp1 - diagVf * conj(Yf * diagV));
    % dSf_dVm = diagVf * conj(Yf * diagVnorm) + conj(diagIf) * temp2;
    % dSt_dVa = 1j * (conj(diagIt) * temp3 - diagVt * conj(Yt * diagV));
    % dSt_dVm = diagVt * conj(Yt * diagVnorm) + conj(diagIt) * temp4;
    %
    % g4_x = [real(dSf_dVa - dSt_dVa); real(dSf_dVm - dSt_dVm)];
    
    %% construct Jacobian of equality (power flow) constraints and transpose it
    % nxyz = nvariable;
    % dg = sparse(2*nb, nxyz);
    % dg(:, [iVa iVm iPg iQg]) = [
    %     real([dSbus_dVa dSbus_dVm]) neg_Cg sparse(nb, ng);  %% P mismatch w.r.t Va, Vm, Pg, Qg
    %     imag([dSbus_dVa dSbus_dVm]) sparse(nb, ng) neg_Cg;  %% Q mismatch w.r.t Va, Vm, Pg, Qg
    %     ];
    % dg = dg';

    %% ====================================================================
    % % form function g first derivatives
    % % Computes partial derivatives of power flows w.r.t. voltage.
    if nl2 > 0
        %% compute partials of Flows w.r.t. V
        if upper(mpopt.opf.flow_lim(1)) == 'I'  %% current
            [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, x0_V_theta);
        else                            %% power
            [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dSbr_dV(branch(il,:), Yf, Yt, x0_V_theta);
        end
        if upper(mpopt.opf.flow_lim(1)) == 'P'  %% real part of flow (active power)
            dFf_dVa = real(dFf_dVa);
            dFf_dVm = real(dFf_dVm);
            dFt_dVa = real(dFt_dVa);
            dFt_dVm = real(dFt_dVm);
            Ff = real(Ff);
            Ft = real(Ft);
        end
        
        %% squared magnitude of flow (of complex power or current, or real power)
        % [df_dVa, df_dVm, dt_dVa, dt_dVm] = ...
        %     dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
        
        %% construct Jacobian of inequality (branch flow) constraints & transpose
        g4_x = [
            dFf_dVa', dFt_dVa';
            dFf_dVm', dFt_dVm'
        ];
        %g4_x = [
        %    df_dVa, df_dVm;
        %    dt_dVa, dt_dVm
        %];
    else
        fprintf('Not supported!');
        success = 0;
        return;
        % g4_x = zeros(2 * nb, nl2);
    end
    
    g_grad_matrix = [eye(ng, ng) zeros(ng, r_num2 - ng);
        zeros(ng, ng), eye(ng, ng), zeros(ng, r_num2 - 2 * ng);
        zeros(2 * nb, ng), zeros(2 * nb, ng), [zeros(nb, nb); eye(nb, nb)], g4_x];
    
    %% ====================================================================
    %% form diag matrix
    L_Z = diag(zLagr./lslack);
    U_W = diag(wLagr./uslack);
    
    %% ====================================================================
    %% form hession matrix
    % form function f first derivatives
    f_grad_matrix = [2 * gencost(:, COST) .* x0_Pg  * 1e4 + gencost(:, COST + 1) * 1e2;
        zeros(ng, 1);
        zeros(2 * nb, 1)];
    
    % form function f second derivatives
    f_hess_matrix = zeros(nvariable, nvariable);
    f_hess_matrix(1:ng, 1:ng) = diag(2 * gencost(:, COST) * 1e4);
    
    %% ====================================================================
    % evaluate Hessian of power balance constraints(h)
    nlam = length(yLagr) / 2;
    lamP = yLagr(1:nlam);
    lamQ = yLagr((1:nlam)+nlam);
    [Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, x0_V_theta, lamP);
    [Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, x0_V_theta, lamQ);
    h_y_hess_matrix = -[
        zeros(2 * ng, nvariable);
        zeros(2 * nb, 2 * ng), real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv])
        ];
    
    %% ====================================================================
    % evaluate Hessian of flow constraints(g)
    nmu = nl2;
    if nmu
        muF = zLagr((1:nmu)+2*ng+nb) + wLagr((1:nmu)+2*ng+nb);
        muT = zLagr((1:nmu)+2*ng+nb+nmu) + wLagr((1:nmu)+2*ng+nb+nmu);
    else    %% keep dimensions of empty matrices/vectors compatible
        muF = zeros(0,1);   %% (required to avoid problems when using Knitro
        muT = zeros(0,1);   %%  on cases with all lines unconstrained)
    end
    % upper(mpopt.opf.flow_lim(1)) != 'I'
    f = branch(il, F_BUS);    %% list of "from" buses
    t = branch(il, T_BUS);    %% list of "to" buses
    Cf = sparse(1:nl2, f, ones(nl2, 1), nl2, nb);     %% connection matrix for line & from buses
    Ct = sparse(1:nl2, t, ones(nl2, 1), nl2, nb);     %% connection matrix for line & to buses
    %[dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, x0_V_theta);
    [Hfaa, Hfav, Hfva, Hfvv] = d2Sbr_dV2(Cf, Yf, x0_V_theta, muF);
    [Htaa, Htav, Htva, Htvv] = d2Sbr_dV2(Ct, Yt, x0_V_theta, muT);
    %     if upper(mpopt.opf.flow_lim(1)) == 'P'    %% real power
    %         [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(real(dSf_dVa), real(dSf_dVm), real(Sf), Cf, Yf, x0_V_theta, muF);
    %         [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(real(dSt_dVa), real(dSt_dVm), real(St), Ct, Yt, x0_V_theta, muT);
    %     else                                      %% apparent power
    %         [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, Yf, x0_V_theta, muF);
    %         [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, Yt, x0_V_theta, muT);
    %     end
    if upper(mpopt.opf.flow_lim(1)) == 'P'  %% real part of flow (active power)
        Hfaa = real(Hfaa);
        Hfav = real(Hfav);
        Hfva = real(Hfva);
        Hfvv = real(Hfvv);
        Htaa = real(Htaa);
        Htav = real(Htav);
        Htva = real(Htva);
        Htvv = real(Htvv);
    end
    g_zw_hess_matrix = [
        zeros(2 * ng, nvariable);
        zeros(2 * nb, 2 * ng), [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv]
        ];
    
    HH_matrix = -f_hess_matrix + h_y_hess_matrix + g_zw_hess_matrix - g_grad_matrix * (L_Z - U_W) * g_grad_matrix';
    
    %% ====================================================================
    %% rebuild Sbus
    Sbus = makeSbus(baseMVA, bus, gen, mpopt, x0_V);  %% net injected power in p.u.
    
    %% evaluate power flow equations
    mis = x0_V_theta .* conj(Ybus * x0_V_theta) - Sbus;
    
    %%----- evaluate constraint function values -----
    %% first, the equality constraints (power flow)
    h_fcn = -[ real(mis);            %% active power mismatch for all buses
        imag(mis) ];          %% reactive power mismatch for all buses
    
    %% then, the inequality constraints (branch flow limits)
    if nl2 > 0
        flow_max = (branch(il, RATE_A)/baseMVA);
        flow_max(flow_max == 0) = Inf;
        if upper(mpopt.opf.flow_lim(1)) == 'I'    %% current magnitude limit, |I|
            fprintf('Not supported!');
            success = 0;
            return;
            % If = Yf * x0_V_theta;
            % It = Yt * x0_V_theta;
            % g_fcn = [ If .* conj(If) - flow_max;    %% branch current limits (from bus)
            %      It .* conj(It) - flow_max ];  %% branch current limits (to bus)
        else
            %% compute branch power flows
            Sf = x0_V_theta(branch(il, F_BUS)) .* conj(Yf * x0_V_theta);  %% complex power injected at "from" bus (p.u.)
            St = x0_V_theta(branch(il, T_BUS)) .* conj(Yt * x0_V_theta);  %% complex power injected at "to" bus (p.u.)
            if upper(mpopt.opf.flow_lim(1)) == 'P'  %% active power limit, P (Pan Wei)
                C_Lz = [x0_Pg - gen(:, PMIN)/baseMVA;
                    x0_Qg - gen(:, QMIN)/baseMVA;
                    x0_V - bus(:, VMIN);
                    real(Sf) - (- flow_max);
                    real(St) - (- flow_max)] - lslack;
                C_Lw = [x0_Pg - gen(:, PMAX)/baseMVA;
                    x0_Qg - gen(:, QMAX)/baseMVA;
                    x0_V - bus(:, VMAX);
                    real(Sf) - flow_max;
                    real(St) - flow_max] + uslack;
            else                                    %% apparent power limit, |S|
                fprintf('Not supported!');
                success = 0;
                return;
                g_fcn = [ Sf .* conj(Sf) - flow_max;      %% branch apparent power limits (from bus)
                    St .* conj(St) - flow_max ];    %% branch apparent power limits (to bus)
            end
        end
    else
        fprintf('Not supported!');
        success = 0;
        return;
        g_fcn = zeros(0,1);
    end
    
    C_Ly = h_fcn;
    C_Ll = diag(lslack)*diag(zLagr)*ones(r_num2, 1) - miu*ones(r_num2, 1);
    C_Lu = diag(uslack)*diag(wLagr)*ones(r_num2, 1) + miu*ones(r_num2, 1);
    C_Lxx = f_grad_matrix - h_grad_matrix * yLagr - g_grad_matrix * (zLagr + wLagr) + g_grad_matrix * (diag(lslack) ^(-1) * (C_Ll + diag(zLagr) * C_Lz) + diag(uslack) ^(-1) * (C_Lu - diag(wLagr) * C_Lw));
    
    
    [deltaEq_xxx, info] = mplinsolve([HH_matrix,h_grad_matrix;h_grad_matrix',zeros(2*nb, 2*nb)], [C_Lxx;-C_Ly], opt.linsolver, []);
    delta_x = deltaEq_xxx(1:nvariable);
    delta_x_x = delta_x';
    delta_y =  deltaEq_xxx(nvariable+1:end);
    % [deltaEq_xxx, info] = mplinsolve(h_grad_matrix, (C_Lxx - HH_matrix * delta_x), opt.linsolver, []);
    
    % [deltaEq_xxx, info] = mplinsolve(eye(r_num2, r_num2), (), opt.linsolver, []);
    delta_u = -C_Lw - g_grad_matrix' * delta_x;
    
    % [deltaEq_xxx, info] = mplinsolve(eye(r_num2, r_num2), (), opt.linsolver, []);
    delta_l = C_Lz + g_grad_matrix' * delta_x;
    
    % [deltaEq_xxx, info] = mplinsolve(eye(r_num2, r_num2), (), opt.linsolver, []);
    delta_z = -diag(lslack)^(-1) * C_Ll - L_Z * delta_l;
    
    % [deltaEq_xxx, info] = mplinsolve(eye(r_num2, r_num2), (), opt.linsolver, []);
    delta_w = -diag(uslack)^(-1) * C_Lu - U_W * delta_u;
    
    %     Eq_AAA = [eye(r_num2, r_num2), L_Z, zeros(r_num2, r_num2), zeros(r_num2, r_num2), zeros(r_num2, nvariable), zeros(r_num2, 2*nb);
    %               zeros(r_num2, r_num2), eye(r_num2, r_num2), zeros(r_num2, r_num2), zeros(r_num2, r_num2), -g_grad_matrix', zeros(r_num2, 2*nb);
    %               zeros(r_num2, r_num2), zeros(r_num2, r_num2), eye(r_num2, r_num2), U_W, zeros(r_num2, nvariable), zeros(r_num2, 2*nb);
    %               zeros(r_num2, r_num2), zeros(r_num2, r_num2), zeros(r_num2, r_num2), eye(r_num2, r_num2), g_grad_matrix', zeros(r_num2, 2*nb);
    %               zeros(nvariable, r_num2), zeros(nvariable, r_num2), zeros(nvariable, r_num2), zeros(nvariable, r_num2), HH_matrix, -h_grad_matrix;
    %               zeros(2*nb, r_num2), zeros(2*nb, r_num2), zeros(2*nb, r_num2), zeros(2*nb, r_num2), -h_grad_matrix', zeros(2*nb, 2*nb);
    %     ];
    %
    %
    %     Eq_BBB = [-diag(lslack)^(-1) * C_Ll;
    %               C_Lz;
    %               -diag(uslack)^(-1) * C_Lu;
    %               -C_Lw;
    %               C_Lxx;
    %               -C_Ly
    %     ];
    %
    %
    %     [deltaEq_xxx, info] = mplinsolve(Eq_AAA, Eq_BBB, opt.linsolver, []);
    %     delta_z = deltaEq_xxx(1:r_num2);
    %     delta_l = deltaEq_xxx(1*r_num2+1 : 2*r_num2);
    %     delta_w = deltaEq_xxx(2*r_num2+1 : 3*r_num2);
    %     delta_u = deltaEq_xxx(3*r_num2+1 : 4*r_num2);
    %     delta_x = deltaEq_xxx(4*r_num2+1: 4*r_num2+nvariable);
    %     delta_y = deltaEq_xxx(4*r_num2+nvariable+1 : end);
    %
    delta_l_below0_idx = find(delta_l < 0);
    delta_u_below0_idx = find(delta_u < 0) ;
    alpha_p = 0.9995 * min([min(-lslack(delta_l_below0_idx) ./ delta_l(delta_l_below0_idx)), min(-uslack(delta_u_below0_idx) ./ delta_u(delta_u_below0_idx)), 1]);
    %     alpha_p = 1;
    %     alpha_d = 1;
    delta_z_below0_idx = find(delta_z < 0);
    delta_w_below0_idx = find(delta_w > 0);
    alpha_d = 0.9995 * min([min(-zLagr(delta_z_below0_idx) ./ delta_z(delta_z_below0_idx)), min(-wLagr(delta_w_below0_idx) ./ delta_w(delta_w_below0_idx)), 1]);
    
    x0 = x0 + alpha_p * delta_x;
    x0(2*ng + 1:2*ng + nb) = x0(2*ng + 1:2*ng + nb) - x0(2*ng + nb);
    lslack = lslack + alpha_p * delta_l;
    uslack = uslack + alpha_p * delta_u;
    yLagr = yLagr + alpha_d * delta_y;
    zLagr = zLagr + alpha_d * delta_z;
    wLagr = wLagr + alpha_d * delta_w;
    
    Gap = lslack' * zLagr - uslack' * wLagr;
    iterCount = iterCount + 1;
end
%
% f_fcn = @(x)opf_costfcn(x, om);
% gh_fcn = @(x)opf_consfcn(x, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
% hess_fcn = @(x, lambda, cost_mult)opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
%
% [x, f, info, Output, Lambda] = ...
%     mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
% success = (info > 0);
%

%%-----  calculate return values  -----
%% update voltages & generator outputs
x0_Pg = x0(1:ng);
x0_Qg = x0(ng + 1:2 * ng);
x0_theta = x0(2 * ng + 1: 2 * ng + nb);
x0_V = x0(2 * ng + nb + 1 : end);
x0_V_theta = x0_V .* exp(1j * x0_theta);
Ibus = Ybus * x0_V_theta;
diagIbus = diag(Ibus);
% form function f
f_value = sum(x0_Pg .^2 .* gencost(:, COST) * 1e4 + x0_Pg .* gencost(:, COST+1) * 1e2 + gencost(:, COST+2));

gen(:, PG) = x0_Pg * baseMVA;
gen(:, QG) = x0_Qg * baseMVA;
gen(:, VG) = x0(2*ng + nb + gen(:, GEN_BUS));
bus(:, VA) = rad2deg(x0_theta);
bus(:, VM) = x0_V;
%% compute branch flows
Sf = x0_V_theta(branch(:, F_BUS)) .* conj(Yf * x0_V_theta);  %% cplx pwr at "from" bus, p.u.
St = x0_V_theta(branch(:, T_BUS)) .* conj(Yt * x0_V_theta);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

%% package up results
nlnN = 2*nb + 2*nl;     %% because muSf and muSt are nl x 1, not nl2 x 1

results = mpc;
% 系统中所有变量，排列顺序为P1~g,Q1~g,theta1~n,V1~n
tmp = x0(2*nb+1:end);
x0(2*nb+1:end) = [];
x0 = [tmp; x0];

[results.bus, results.branch, results.gen, results.x, results.f] = deal(bus, branch, gen, x0, f_value);
success = 1;
raw = struct('xr', x0, 'info', info);

plot(1:size(Gap_arr,2),Gap_arr,'-*');
set(gca,'YScale','log');
xlabel('迭代次数'),ylabel('Gap');
title('最优潮流内点法收敛特性');