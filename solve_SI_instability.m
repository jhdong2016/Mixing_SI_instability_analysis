function R = solve_SI_instability(p, opt)
%SOLVE_SI_INSTABILITY Discretize and solve the manuscript's linear problem.
%   R = solve_SI_instability(p,opt)
%
%   This is an eigenvalue SOLVER, not a fit formula or a stored-data reader.
%   The companion main_SI_instability.m contains the editable parameters.
%   All helpers are local to this file; no third-party .m files are required.
%
%   REQUIRED PHYSICAL PARAMETERS (one scalar case per call)
%     p.Ri       gradient Richardson number, > 0
%     p.Ev       vertical Ekman number, >= 0
%     p.delta    H/L = f/S0, > 0
%     p.Pr       turbulent Prandtl number, > 0; Eb = Ev/Pr
%
%   OPTIONAL MODEL PARAMETERS
%     p.buoyancy_bc       'neumann' (default) or 'dirichlet'
%     p.horizontal_ratio  nu_h/nu_v = kappa_h/kappa_v (default 0)
%                         This is a DIMENSIONAL diffusivity ratio. In the
%                         manuscript units the extra term is -Ev*r*delta^2*K^2.
%     p.U_offset          U(z)=z+U_offset; default 0.5. Does not change U_z.
%
%   NUMERICAL OPTIONS (all settable in the companion main program)
%     opt.task            'si' | 'point' | 'spectrum'
%     opt.method          'chebyshev' (default) | 'staggered'
%     opt.N               96 Chebyshev intervals / 80 staggered cells
%     opt.scan_N          optional coarse-scan resolution; [] uses opt.N
%                         local refinement/final eigenfunctions use opt.N
%     opt.K_grid          positive total horizontal wavenumbers
%     opt.phi_deg         sampled angles for 'spectrum'; k=K*sin(phi),
%                         ell=K*cos(phi). 'si' always uses exactly k=0.
%     opt.k, opt.ell      inputs for 'point' (manuscript convention)
%     opt.refine          locally refine K at sampled angles, default true
%     opt.logK_tolerance  fminbnd tolerance in log(K), default 1e-6
%     opt.growth_tol      detection threshold for Re(sigma), default 1e-6
%     opt.qz_tol          homogeneous denominator filter, default 1e-9
%     opt.check_N         optional same-wavenumber resolution check, [] off
%     opt.check_peak      also reoptimize the peak bracket at check_N
%     opt.check_no_growth confirm marginal/non-growing Cheb scans with FD
%     opt.confirm_N      staggered cells for that check, default 160
%     opt.confirm_trigger 1e-5 by default; see README for interpretation
%     opt.return_mode     return the leading eigenfunction, default true
%     opt.return_matrices return the final A,B matrices, default false
%
%   MODEL AND CONVENTIONS
%     z in [-1,0], U=z+1/2, V=W=0, B_z=1, B_y=-1/Ri.
%     Perturbations: exp(sigma*t + 1i*(k*x + ell*y)); Q=sigma+1i*k*U.
%       Qu+w-v = -Ri*1i*k*p + Ev*D^2*u
%       Qv+u   = -Ri*1i*ell*p + Ev*D^2*v
%       delta^2*Qw = -Ri*Dp + Ri*b + Ev*delta^2*D^2*w
%       Qb-v/Ri+w = (Ev/Pr)*D^2*b
%       1i*k*u+1i*ell*v+Dw = 0
%     Finite-mixing BCs: w=0, Du=Dv=Db=0 (unless buoyancy BC changed).
%     These BCs act on perturbations, not on the maintained basic state.
%     sigma and its real part are NONDIMENSIONAL; do not divide by f again.
%
%   METHODS
%     'chebyshev' ports the archived Figure-3 generalized (zeta,w,beta)
%       collocation matrix, including its internal coordinate rotation and
%       energy scaling. QZ excludes infinite boundary-constraint eigenvalues.
%     'staggered' ports the primitive, energy-compatible pressure projection
%       used for Figures 1, 2, and 4. Pressure is recovered by least squares.
%     At Ev=0 the requested Cheb method switches EXPLICITLY in R.notes to the
%       staggered method: viscous stress/scalar BCs must not be imposed on
%       the reduced-order inviscid problem. No inviscid strict-SI selected
%       wavelength is reported from a finite scan.
%
%   SEARCH SCOPE
%     'si': scans ell=K>0 at k=0, then refines positive local maxima in log K.
%     'spectrum': scans (K,phi) and refines K only at the sampled phi values.
%                 The angular maximum remains a SAMPLED-ANGLE maximum.
%     'point': returns all finite eigenvalues at the requested (k,ell).
%     No theoretically predicted wavelength is used to locate a peak.
%
%   OUTPUTS
%     R.peak: sigma,growth,k,ell,K,phi_deg,lambda_over_L,status,has_growth,
%       selected_wavelength_defined. For unresolved/edge/no-growth maxima,
%       lambda_over_L is NaN; candidate_lambda_over_L keeps the tested scale.
%     R.scan: K,phi_deg,k,ell,sigma_lead,growth (including neutral/decaying).
%     R.eigenvalues / unstable_eigenvalues: all finite / growing eigenvalues
%       at the final candidate point; only the leading eigenfunction is saved.
%     R.mode: u,v,w,b,p,z plus residual/energy diagnostics when requested.
%     R.check: separate fixed-wavenumber / optimized-peak resolution results.
%     R.confirmation: separate primitive-solver check near zero growth.
%
%   LIMITATION: a small algebraic residual is not proof of spatial convergence.
%   Increasing N and expanding the scan are needed when changing parameters.

if nargin < 1 || ~isstruct(p), error('SI:Input','p must be a structure.'); end
if nargin < 2, opt = struct(); end
[p,opt] = validate_inputs(p,opt);
R = struct('parameters',p,'options',opt,'method_requested',opt.method, ...
    'method_used',opt.method,'notes',{{}},'confirmation',struct(), ...
    'check',struct(),'mode',struct(),'eigenvalues',[],'unstable_eigenvalues',[]);
if p.Ev == 0 && strcmp(opt.method,'chebyshev')
    opt.method = 'staggered'; R.method_used = opt.method;
    R.notes{end+1} = ['Ev=0: used staggered primitive equations, not the ', ...
        'finite-viscosity Chebyshev boundary rows.'];
end
O = make_operators(opt.N,opt.method,p.buoyancy_bc);

if strcmp(opt.task,'point')
    kk=opt.k; ll=opt.ell; KK=hypot(kk,ll); angle=atan2(kk,ll)*180/pi;
    scanSigma=mode_value(p,opt,O,kk,ll);
    R.scan=struct('K',KK,'phi_deg',angle,'k',kk,'ell',ll, ...
        'sigma_lead',scanSigma,'growth',real(scanSigma));
    candidate=struct('K',KK,'phi',angle,'sigma',scanSigma,'bracket',[], ...
        'radial_boundary',false,'optimizer_exitflag',NaN);
else
    K=opt.K_grid(:); phi=opt.phi_deg(:).';
    if strcmp(opt.task,'si'), phi=0; end
    ns=numel(K); na=numel(phi);
    scanopt=opt;
    if ~isempty(opt.scan_N), scanopt.N=opt.scan_N; end
    Oscan=make_operators(scanopt.N,scanopt.method,p.buoyancy_bc);
    S=complex(zeros(ns,na)); kmap=zeros(ns,na); lmap=kmap;
    candidates=struct('K',{},'phi',{},'sigma',{},'bracket',{}, ...
        'radial_boundary',{},'optimizer_exitflag',{});
    fmopt=optimset('Display','off','TolX',opt.logK_tolerance, ...
        'MaxFunEvals',opt.max_refine_evaluations,'MaxIter',opt.max_refine_evaluations);
    for ia=1:na
        [kc,lc]=wavevector(1,phi(ia));
        for ik=1:ns
            kmap(ik,ia)=K(ik)*kc; lmap(ik,ia)=K(ik)*lc;
            S(ik,ia)=mode_value(p,scanopt,Oscan,kmap(ik,ia),lmap(ik,ia));
        end
        g=real(S(:,ia));
        % Retain every scanned point; do not discard stable/neutral eigenvalues.
        [~,ig]=max(g);
        for jq=unique([1,ig,ns])
            sh=mode_value(p,opt,O,K(jq)*kc,K(jq)*lc);
            candidates(end+1)=make_candidate(K(jq),phi(ia),sh,[], ...
                jq==1 || jq==ns,NaN); %#ok<AGROW>
        end
        if opt.refine
            for ik=2:ns-1
                if g(ik)>opt.growth_tol && g(ik)>=g(ik-1) && g(ik)>=g(ik+1)
                    br=log(K([ik-1,ik+1])).';
                    objective=@(x) -real(mode_value(p,opt,O,exp(x)*kc,exp(x)*lc));
                    [xx,~,exitflag]=fminbnd(objective,br(1),br(2),fmopt);
                    sr=mode_value(p,opt,O,exp(xx)*kc,exp(xx)*lc);
                    candidates(end+1)=make_candidate(exp(xx),phi(ia),sr,br, ...
                        false,exitflag); %#ok<AGROW>
                end
            end
        end
    end
    [~,best]=max(real([candidates.sigma])); candidate=candidates(best);
    R.scan=struct('K',K,'phi_deg',phi,'k',kmap,'ell',lmap, ...
        'sigma_lead',S,'growth',real(S),'N',scanopt.N);
    R.local_candidates=candidates;
end

if strcmp(opt.task,'point')
    kstar=opt.k; lstar=opt.ell;
else
    [kstar,lstar]=wavevector(candidate.K,candidate.phi);
end
P=solve_point(p,opt,O,kstar,lstar,opt.return_mode,opt.return_matrices);
R.eigenvalues=P.eigenvalues;
R.unstable_eigenvalues=P.eigenvalues(real(P.eigenvalues)>opt.growth_tol);
R.n_finite=numel(P.eigenvalues); R.n_unstable=numel(R.unstable_eigenvalues);
R.n_qz_excluded=P.n_qz_excluded;
if opt.return_mode, R.mode=P.mode; end
if opt.return_matrices, R.A=P.A; R.B=P.B; end
sigma=P.sigma; growing=real(sigma)>opt.growth_tol;
Kboundary=candidate.radial_boundary;
angleEdge=false;
if strcmp(opt.task,'spectrum') && numel(opt.phi_deg)>1
    angleEdge=(candidate.phi==min(opt.phi_deg) || candidate.phi==max(opt.phi_deg));
end

% Marginal Chebyshev eigenvalues can have tiny roundoff-induced positive real
% parts. Do not silently interpret them as growing modes or as exact zeros.
ambiguous=false;
if strcmp(opt.method,'chebyshev') && opt.check_no_growth && real(sigma)<=opt.confirm_trigger
    of=opt; of.method='staggered'; of.N=opt.confirm_N;
    Of=make_operators(of.N,of.method,p.buoyancy_bc);
    cg=nan(size(R.scan.k));
    for q=1:numel(cg)
        cg(q)=real(mode_value(p,of,Of,R.scan.k(q),R.scan.ell(q)));
    end
    [gm,ii]=max(cg(:));
    g_at_candidate=real(mode_value(p,of,Of,kstar,lstar));
    gm=max(gm,g_at_candidate);
    R.confirmation=struct('method','staggered','N',of.N,'scan_growth',cg, ...
        'scan_max_growth',gm,'growth_at_spectral_candidate',g_at_candidate, ...
        'k_at_scan_max',R.scan.k(ii),'ell_at_scan_max',R.scan.ell(ii), ...
        'agrees_on_growth',growing==(gm>opt.growth_tol));
    ambiguous=~R.confirmation.agrees_on_growth;
    if ambiguous
        R.notes{end+1}=['Near-neutral classification differs between methods. ', ...
            'Increase N and inspect the spectrum; no selected wavelength is certified.'];
    end
end

optimizer_failed=~isnan(candidate.optimizer_exitflag) && candidate.optimizer_exitflag<=0;
local_edge=false;
if numel(candidate.bracket)==2
    local_edge=min(log(candidate.K)-candidate.bracket(1),candidate.bracket(2)-log(candidate.K)) ...
        <=10*opt.logK_tolerance;
end
selected=growing && ~Kboundary && ~ambiguous && ~optimizer_failed && ~local_edge && ~strcmp(opt.task,'point');
if strcmp(opt.task,'si') && p.Ev==0
    selected=false;
    R.notes{end+1}=['Inviscid strict SI: a finite scan cannot identify a physical ', ...
        'fastest-growing wavelength. Use the scanned growth as a finite-K diagnostic.'];
end
if strcmp(opt.task,'point')
    if growing, status='Growing mode at specified wavenumber';
    else, status='No detected growth at specified wavenumber'; end
elseif ambiguous
    status='Unresolved near-neutral classification; compare discretizations';
elseif ~growing
    status='No detected positive growth in the searched range';
elseif strcmp(opt.task,'si') && p.Ev==0
    status='Inviscid strict SI; no finite selected wavelength assigned';
elseif optimizer_failed || local_edge
    status='Peak refinement unresolved; refine or expand the numerical bracket';
elseif Kboundary
    status='Growing maximum on K scan boundary; expand the range';
elseif strcmp(opt.task,'spectrum')
    status='Growing interior-K maximum over sampled angles';
else
    status='Growing interior strict-SI maximum';
end
lambdaSelected=NaN;
if selected, lambdaSelected=2*pi/candidate.K; end
R.peak=struct('sigma',sigma,'growth',real(sigma),'k',kstar,'ell',lstar, ...
    'K',candidate.K,'phi_deg',candidate.phi,'lambda_over_L',lambdaSelected, ...
    'candidate_lambda_over_L',2*pi/candidate.K,'has_growth',growing, ...
    'selected_wavelength_defined',selected,'radial_boundary',Kboundary, ...
    'angle_at_scan_edge',angleEdge,'local_bracket_edge',local_edge,'status',status, ...
    'logK_bracket',candidate.bracket,'optimizer_exitflag',candidate.optimizer_exitflag);
if strcmp(opt.task,'spectrum')
    R.notes{end+1}=['K is locally optimized at sampled phi only. Refine phi_deg ', ...
        'to test an oblique maximum; this is not continuous 2-D optimization.'];
end

% A fixed-wavenumber check and a reoptimized-wave-number check are different.
if ~isempty(opt.check_N)
    oc=opt; oc.N=opt.check_N;
    Oc=make_operators(oc.N,oc.method,p.buoyancy_bc);
    sc=mode_value(p,oc,Oc,kstar,lstar);
    R.check=struct('method',oc.method,'N',oc.N,'same_K_sigma',sc, ...
        'same_K_growth_difference',real(sc)-real(sigma), ...
        'same_K_relative_growth_difference_pct',NaN, ...
        'reoptimized',false,'reoptimized_K',NaN,'reoptimized_sigma',NaN, ...
        'relative_wavelength_difference_pct',NaN);
    if abs(real(sigma))>opt.growth_tol
        R.check.same_K_relative_growth_difference_pct=100*(real(sc)/real(sigma)-1);
    end
    if opt.check_peak && selected && numel(candidate.bracket)==2
        br=candidate.bracket;
        [kd,ld]=wavevector(1,candidate.phi);
        fun=@(x) -real(mode_value(p,oc,Oc,exp(x)*kd,exp(x)*ld));
        fmopt=optimset('Display','off','TolX',opt.logK_tolerance, ...
            'MaxFunEvals',opt.max_refine_evaluations,'MaxIter',opt.max_refine_evaluations);
        [xc,~,ef]=fminbnd(fun,br(1),br(2),fmopt);
        sr=mode_value(p,oc,Oc,exp(xc)*kd,exp(xc)*ld);
        R.check.reoptimized=true; R.check.optimizer_exitflag=ef;
        R.check.reoptimized_K=exp(xc); R.check.reoptimized_sigma=sr;
        R.check.relative_wavelength_difference_pct=100*(candidate.K/exp(xc)-1);
    end
end
end

function c=make_candidate(K,phi,sigma,br,edge,ef)
c=struct('K',K,'phi',phi,'sigma',sigma,'bracket',br, ...
    'radial_boundary',edge,'optimizer_exitflag',ef);
end

function [k,ell]=wavevector(K,phi)
k=K*sind(phi); ell=K*cosd(phi);
% Ensure the named axes are exact, rather than tiny nonzero trigonometric values.
if abs(sind(phi))<1e-14, k=0; end
if abs(cosd(phi))<1e-14, ell=0; end
end

function s=mode_value(p,opt,O,k,ell)
P=solve_point(p,opt,O,k,ell,false,false); s=P.sigma;
end

function P=solve_point(p,opt,O,kman,ellman,needMode,needMatrices)
% Archive working frame: x_c=y, y_c=-x; all lengths scaled by H.
% u_man=-delta*v_c, v_man=delta*u_c, w_man=w_c, b_man=beta/sqrt(S2).
% k_c=delta*ell_man, ell_c=-delta*k_man; no change of physical orientation.
kc=p.delta*ellman; lc=-p.delta*kman; K2=kc^2+lc^2;
if K2<=0, error('SI:ZeroWave','k=ell=0 is excluded.'); end
M=-1/p.delta; S2=p.Ri/p.delta^2; rootS=sqrt(S2);
nu=p.Ev; kap=nu/p.Pr; rh=p.horizontal_ratio;
if strcmp(opt.method,'staggered')
    N=opt.N; I=O.I; J=O.J; D=O.D; T=O.T;
    Vc=M*(O.zc+p.U_offset); Vf=M*(O.zf+p.U_offset);
    Hc=nu*(O.Lc-rh*K2*I)-1i*lc*diag(Vc);
    Hf=nu*(O.Lf-rh*K2*J)-1i*lc*diag(Vf);
    Hb=kap*(O.Lb-rh*K2*I)-1i*lc*diag(Vc);
    nv=3*N-1; nr=2*N-1; iu=1:N; iv=N+(1:N); iw=2*N+(1:N-1);
    F=complex(zeros(nv)); F(iu,iu)=Hc; F(iu,iv)=I;
    F(iv,iu)=-I; F(iv,iv)=Hc; F(iv,iw)=-M*T; F(iw,iw)=Hf;
    Bb=complex(zeros(nv,N)); Bb(iw,:)=rootS*T';
    Cb=complex(zeros(N,nv)); Cb(:,iu)=-M/rootS*I; Cb(:,iw)=-rootS*T;
    Kc=sqrt(K2); co=kc/Kc; si=lc/Kc;
    % Compatible orthonormal null-space basis (no explicit matrix inverse).
    Ch=chol(J+(D'*D)/K2); Tf=Ch\J; DT=D*Tf;
    Z=complex(zeros(nv,nr)); Z(iu,1:N)=1i*si*I; Z(iv,1:N)=-1i*co*I;
    Z(iu,N+1:end)=1i*co/Kc*DT; Z(iv,N+1:end)=1i*si/Kc*DT;
    Z(iw,N+1:end)=Tf;
    A=[Z'*F*Z, Z'*Bb; Cb*Z, Hb]; B=eye(size(A));
    if needMode
        [EV,ED]=eig(A); ev=diag(ED);
    else
        ev=eig(A); EV=[];
    end
    if any(~isfinite(ev)), error('SI:FiniteSpectrum','Nonfinite projected eigenvalue.'); end
    P.n_qz_excluded=0;
    [~,ix]=max(real(ev)); s=ev(ix);
    if needMode
        x=EV(:,ix); vel=Z*x(1:nr); beta=x(nr+1:end);
        div=[1i*kc*I,1i*lc*I,D]; grad=-div';
        force=F*vel+Bb*beta;
        pc=grad\(force-s*vel);  % Rectangular least squares, not normal equations.
        rm=s*vel-force+grad*pc; rb=s*beta-Cb*vel-Hb*beta;
        scale=norm(s*[vel;beta])+norm([force;Cb*vel+Hb*beta])+norm(grad*pc)+realmin;
        rr=norm([rm;rb])/scale;
        uc=vel(iu); vc=vel(iv); wc=vel(iw);
        E=.5*real(vel'*vel+beta'*beta);
        PS=-M*real(vc'*(T*wc)); PF=-M/rootS*real(beta'*uc);
        DV=-nu*real(uc'*(O.Lc-rh*K2*I)*uc+vc'*(O.Lc-rh*K2*I)*vc+wc'*(O.Lf-rh*K2*J)*wc);
        DB=-kap*real(beta'*(O.Lb-rh*K2*I)*beta);
        mz2=-real(uc'*O.Lc*uc+vc'*O.Lc*vc+wc'*O.Lf*wc)/max(real(vel'*vel),realmin);
        enres=abs(2*real(s)*E-PS-PF+DV+DB)/(abs(2*real(s)*E)+abs(PS)+abs(PF)+abs(DV)+abs(DB)+realmin);
        mode=struct('z',O.zc,'z_w',[-1;O.zf;0], ...
            'u',-p.delta*vc,'v',p.delta*uc,'w',[0;wc;0], ...
            'b',beta/rootS,'p',pc/S2,'zeta',1i*(kman*p.delta*uc+ellman*p.delta*vc), ...
            'primitive_residual',rr,'continuity_residual',norm(div*vel)/(norm(div,'fro')*norm(vel)+realmin), ...
            'energy_residual',enres,'energy_working_units',E, ...
            'kinetic_energy_working_units',.5*real(vel'*vel), ...
            'potential_energy_working_units',.5*real(beta'*beta), ...
            'shear_rate',PS/(2*E),'frontal_rate',PF/(2*E), ...
            'viscous_loss_rate',DV/(2*E),'buoyancy_loss_rate',DB/(2*E), ...
            'm_eff2',mz2,'energy_note','Common uniform cell width is omitted; ratios are unaffected.');
        mode.eigen_residual=norm(A*x-s*x)/((norm(A,'fro')+abs(s))*norm(x)+realmin);
        mode.discrete_eigenvector=x;
        mode.discrete_coordinates='[a; beta], vel_working=Z*a; projected pressure-free system';
        P.mode=normalize_mode(mode);
    end
else
    n=opt.N+1; I=O.I; D=O.D; D2=O.D2;
    Delta=D2-K2*I; Diff=D2-rh*K2*I; V=diag(M*(O.z+p.U_offset));
    iz=1:n; iw=n+(1:n); ib=2*n+(1:n);
    A=complex(zeros(3*n)); B=complex(zeros(3*n));
    A(iz,iz)=nu*Diff-1i*lc*V;
    A(iz,iw)=D-1i*kc*M*I;
    A(iw,iz)=-D;
    A(iw,iw)=nu*(Diff*Delta)-1i*lc*(V*Delta);
    A(iw,ib)=-K2*rootS*I;
    A(ib,iz)=-1i*M*lc/(K2*rootS)*I;
    A(ib,iw)=-rootS*I-1i*M*kc/(K2*rootS)*D;
    A(ib,ib)=kap*Diff-1i*lc*V;
    B(iz,iz)=I; B(iw,iw)=Delta; B(ib,ib)=I;
    % Boundary-row replacement: D*zeta=0, w=0, D^2*w=0, Db=0.
    for j=[1,n]
        rows=[iz(j),iw(j),ib(j)]; A(rows,:)=0; B(rows,:)=0;
        A(iz(j),iz)=D(j,:); A(iw(j),iw(j))=1;
        if strcmp(p.buoyancy_bc,'neumann'), A(ib(j),ib)=D(j,:);
        else, A(ib(j),ib(j))=1; end
    end
    rows=[iw(2),iw(n-1)]; A(rows,:)=0; B(rows,:)=0;
    A(iw(2),iw)=D2(1,:); A(iw(n-1),iw)=D2(n,:);
    % B is singular because the boundary constraints have no time derivative.
    % Never use inv(B)*A or B\A. QZ supplies the homogeneous eigenvalue pair.
    if needMode
        [AA,BB,~,~,EV]=qz(A,B,'complex');
    else
        [AA,BB]=qz(A,B,'complex'); EV=[];
    end
    aa=diag(AA); bb=diag(BB);
    keep=abs(bb)>opt.qz_tol*(abs(aa)+abs(bb));
    keep=keep & isfinite(aa) & isfinite(bb);
    P.n_qz_excluded=sum(~keep);
    ev=aa(keep)./bb(keep);
    if isempty(ev), error('SI:NoFiniteSpectrum','No finite eigenvalues passed the QZ filter.'); end
    [~,ix]=max(real(ev)); s=ev(ix);
    if needMode
        raw=find(keep); x=EV(:,raw(ix));
        zz=x(iz); ww=x(iw); beta=x(ib);
        % k_c*u_c+ell_c*v_c=i*Dw, zeta=i*k_c*v_c-i*ell_c*u_c.
        uc=1i*(kc*(D*ww)+lc*zz)/K2;
        vc=1i*(lc*(D*ww)-kc*zz)/K2;
        H=nu*Diff-1i*lc*V;
        force=[H*uc+vc-s*uc; H*vc-uc-M*ww-s*vc; H*ww+rootS*beta-s*ww];
        grad=[1i*kc*I;1i*lc*I;D];
        % Recover pressure from ALL collocation equations by least squares.
        % The resulting primitive residual is diagnostic, not a mode filter.
        pc=grad\force;
        rm=grad*pc-force;
        rb=s*beta-(kap*Diff-1i*lc*V)*beta+M/rootS*uc+rootS*ww;
        rr=norm([rm;rb])/(norm(force)+norm(grad*pc)+norm(s*beta)+ ...
            norm((kap*Diff-1i*lc*V)*beta)+norm(M/rootS*uc+rootS*ww)+realmin);
        wt=O.weights;
        sq=@(v) real(sum(wt.*abs(v).^2));
        prod=@(u,v) real(sum(wt.*conj(u).*v));
        EK=.5*(sq(uc)+sq(vc)+sq(ww)); EP=.5*sq(beta); E=EK+EP;
        PS=-M*prod(vc,ww); PF=-M/rootS*prod(beta,uc);
        vert=sq(D*uc)+sq(D*vc)+sq(D*ww);
        DV=nu*(vert+rh*K2*(2*EK)); DB=kap*(sq(D*beta)+rh*K2*(2*EP));
        enres=abs(2*real(s)*E-PS-PF+DV+DB)/(abs(2*real(s)*E)+abs(PS)+abs(PF)+abs(DV)+abs(DB)+realmin);
        mode=struct('z',O.z,'z_w',O.z,'u',-p.delta*vc,'v',p.delta*uc, ...
            'w',ww,'b',beta/rootS,'p',pc/S2,'zeta',zz, ...
            'primitive_residual',rr, ...
            'continuity_residual',norm(1i*kc*uc+1i*lc*vc+D*ww)/( ...
                abs(kc)*norm(uc)+abs(lc)*norm(vc)+norm(D,'fro')*norm(ww)+realmin), ...
            'energy_residual',enres,'energy_working_units',E, ...
            'kinetic_energy_working_units',EK,'potential_energy_working_units',EP, ...
            'shear_rate',PS/(2*E),'frontal_rate',PF/(2*E), ...
            'viscous_loss_rate',DV/(2*E),'buoyancy_loss_rate',DB/(2*E), ...
            'm_eff2',vert/max(2*EK,realmin), ...
            'energy_note','Clenshaw-Curtis quadrature; collocation energy closure is diagnostic.');
        mode.eigen_residual=norm(A*x-s*(B*x))/((norm(A,'fro')+abs(s)*norm(B,'fro'))*norm(x)+realmin);
        mode.boundary_residual=norm(A([iz([1,n]),iw([1,2,n-1,n]),ib([1,n])],:)*x)/(norm(x)+realmin);
        mode.discrete_eigenvector=x;
        mode.discrete_coordinates='[zeta; w; beta] in the documented working frame';
        P.mode=normalize_mode(mode);
    end
end
[~,ord]=sort(real(ev),'descend');
P.sigma=s; P.eigenvalues=ev(ord);
if needMatrices, P.A=A; P.B=B; end
end

function m=normalize_mode(m)
% Put all primitive amplitudes on the SAME arbitrary complex scale. Diagnostics
% above are amplitude independent. Retain the applied scale and raw norm.
v=[m.u;m.v;m.w]; [amp,j]=max(abs(v));
if amp>0
    factor=exp(-1i*angle(v(j)))/amp;
    for name={'u','v','w','b','p','zeta','discrete_eigenvector'}
        m.(name{1})=m.(name{1})*factor;
    end
    for name={'energy_working_units','kinetic_energy_working_units','potential_energy_working_units'}
        m.(name{1})=m.(name{1})*abs(factor)^2;
    end
    m.amplitude_multiplier=factor;
else
    m.amplitude_multiplier=1;
end
m.normalization='Largest primitive-velocity component has unit real amplitude.';
end

function O=make_operators(N,method,bc)
if strcmp(method,'staggered')
    h=1/N; D=zeros(N,N-1); T=D;
    for j=1:N-1
        D(j,j)=1/h; D(j+1,j)=-1/h; T(j,j)=.5; T(j+1,j)=.5;
    end
    Lc=-D*D'; Lf=-D'*D; Lb=Lc;
    if strcmp(bc,'dirichlet'), Lb(1,1)=Lb(1,1)-2/h^2; Lb(end,end)=Lb(end,end)-2/h^2; end
    O=struct('D',D,'T',T,'Lc',Lc,'Lf',Lf,'Lb',Lb,'I',eye(N),'J',eye(N-1), ...
        'zc',-1+((1:N)'-.5)*h,'zf',-1+(1:N-1)'*h,'h',h);
else
    j=(0:N)'; xx=cos(pi*j/N); cc=[2;ones(N-1,1);2].*(-1).^j;
    dx=xx-xx'; D=(cc*(1./cc)')./(dx+eye(N+1));
    D=D-diag(sum(D,2)); D=2*D;
    O=struct('z',(xx-1)/2,'D',D,'D2',D*D,'I',eye(N+1), ...
        'weights',clenshaw_curtis_weights(N)/2);
end
end

function w=clenshaw_curtis_weights(N)
% Quadrature weights on [-1,1], with the same cos(pi*j/N) ordering.
theta=pi*(0:N)'/N; w=zeros(N+1,1); ii=2:N; v=ones(N-1,1);
if mod(N,2)==0
    w(1)=1/(N^2-1); w(end)=w(1);
    for k=1:N/2-1, v=v-2*cos(2*k*theta(ii))/(4*k^2-1); end
    v=v-cos(N*theta(ii))/(N^2-1);
else
    w(1)=1/N^2; w(end)=w(1);
    for k=1:(N-1)/2, v=v-2*cos(2*k*theta(ii))/(4*k^2-1); end
end
w(ii)=2*v/N;
end

function [p,o]=validate_inputs(p,o)
for name={'Ri','Ev','delta','Pr'}
    if ~isfield(p,name{1}), error('SI:Parameter','Missing p.%s.',name{1}); end
    a=p.(name{1});
    if ~(isnumeric(a)&&isreal(a)&&isscalar(a)&&isfinite(a))
        error('SI:Parameter','p.%s must be a finite real scalar.',name{1});
    end
end
if p.Ri<=0 || p.delta<=0 || p.Pr<=0 || p.Ev<0
    error('SI:Parameter','Require Ri>0, delta>0, Pr>0, Ev>=0.');
end
p=fill(p,'buoyancy_bc','neumann'); p=fill(p,'horizontal_ratio',0); p=fill(p,'U_offset',.5);
p.buoyancy_bc=lower(char(p.buoyancy_bc));
if ~any(strcmp(p.buoyancy_bc,{'neumann','dirichlet'})), error('SI:BC','Unknown buoyancy BC.'); end
validateattributes(p.horizontal_ratio,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(p.U_offset,{'numeric'},{'real','finite','scalar'});
if ~isstruct(o), error('SI:Input','opt must be a structure.'); end
o=fill(o,'task','si'); o=fill(o,'method','chebyshev');
o.task=lower(char(o.task)); o.method=lower(char(o.method));
if ~any(strcmp(o.task,{'si','point','spectrum'})), error('SI:Task','task must be si, point, or spectrum.'); end
if ~any(strcmp(o.method,{'chebyshev','staggered'})), error('SI:Method','Unknown discretization.'); end
if strcmp(o.method,'staggered'), defaultN=80; else, defaultN=96; end
o=fill(o,'N',defaultN); o=fill(o,'scan_N',[]); o=fill(o,'K_grid',logspace(-1,log10(350),41));
o=fill(o,'phi_deg',0); o=fill(o,'k',0); o=fill(o,'ell',25);
o=fill(o,'refine',true); o=fill(o,'logK_tolerance',1e-6);
o=fill(o,'growth_tol',1e-6); o=fill(o,'qz_tol',1e-9);
o=fill(o,'max_refine_evaluations',120); o=fill(o,'check_N',[]);
o=fill(o,'check_peak',true); o=fill(o,'check_no_growth',true);
o=fill(o,'confirm_N',160); o=fill(o,'confirm_trigger',1e-5);
o=fill(o,'return_mode',true); o=fill(o,'return_matrices',false);
validateattributes(o.N,{'numeric'},{'real','finite','scalar','integer','>=',8});
validateattributes(o.confirm_N,{'numeric'},{'real','finite','scalar','integer','>=',8});
for name={'logK_tolerance','growth_tol','qz_tol','confirm_trigger'}
    validateattributes(o.(name{1}),{'numeric'},{'real','finite','scalar','positive'});
end
if ~isempty(o.scan_N), validateattributes(o.scan_N,{'numeric'},{'real','finite','scalar','integer','>=',8}); end
if ~isempty(o.check_N), validateattributes(o.check_N,{'numeric'},{'real','finite','scalar','integer','>=',8}); end
validateattributes(o.max_refine_evaluations,{'numeric'},{'real','finite','scalar','integer','>=',10});
for name={'refine','check_peak','check_no_growth','return_mode','return_matrices'}
    val=o.(name{1});
    if ~(isscalar(val)&&(islogical(val)||(isnumeric(val)&&any(val==[0,1]))))
        error('SI:Option','opt.%s must be true or false.',name{1});
    end
    o.(name{1})=logical(val);
end
if strcmp(o.task,'point')
    validateattributes(o.k,{'numeric'},{'real','finite','scalar'});
    validateattributes(o.ell,{'numeric'},{'real','finite','scalar'});
    if hypot(o.k,o.ell)==0, error('SI:ZeroWave','k=ell=0 is excluded.'); end
else
    validateattributes(o.K_grid,{'numeric'},{'real','finite','vector','positive'});
    o.K_grid=sort(unique(o.K_grid(:)));
    if numel(o.K_grid)<3, error('SI:Scan','At least 3 positive K values are required.'); end
    validateattributes(o.phi_deg,{'numeric'},{'real','finite','vector','nonempty'});
    o.phi_deg=sort(unique(o.phi_deg(:).'));
    if strcmp(o.task,'si'), o.phi_deg=0; end
end
end

function s=fill(s,name,value)
if ~isfield(s,name), s.(name)=value; end
end
