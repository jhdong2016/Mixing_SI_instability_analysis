%% MAIN_SI_INSTABILITY -- editable parameters and execution
% Put this file and solve_SI_instability.m in the same folder, then run:
%     main_SI_instability
% No archived MAT/CSV data are required to solve the equations.
% All dimensional conversions are intentionally outside the solver.
clear; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here=pwd; end
addpath(here);

%% 1. Physical nondimensional parameters -- EDIT HERE
Ri_values = 0.70;
Ev_values = [1e-4, 1e-3, 3e-3];
p.delta = 0.10;                         % delta=H/L=f/S0; strictly >0
p.Pr = 1.00;                            % Eb=Ev/Pr
p.buoyancy_bc = 'neumann';               % 'neumann': Db=0; 'dirichlet': b=0
p.horizontal_ratio = 0;                 % 0: vertical-only; 1: equal dimensional h/v coefficients
p.U_offset = 0.5;                       % prescribed U(z)=z+0.5, NOT a TTW adjustment

% To recalculate the full current Figure 3 sweep, replace the two vectors above:
% Ri_values = [0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80];
% Ev_values = [1e-5,3e-5,1e-4,3e-4,1e-3,1.5e-3,2e-3,3e-3,4e-3,6e-3,1e-2];

%% 2. Discretization and search -- EDIT HERE
opt.task = 'si';                        % 'si' | 'point' | 'spectrum'
opt.method = 'staggered';               % all main-text numerical results, including Figure 3
opt.N = 80;                            % number of staggered-grid cells
opt.scan_N = [];                       % use N=80 for both the scan and final peak refinement
% Independent Chebyshev check only (not the main-text data source):
% opt.method = 'chebyshev'; opt.N = 96; opt.scan_N = 64;

opt.K_grid = logspace(-1,log10(350),41); % K=sqrt(k^2+ell^2), positive
opt.phi_deg = -90:5:90;                 % only used when task='spectrum'
% For a more focused near-SI test, e.g.: opt.phi_deg = -10:0.5:10;
% In task='si', k=0 exactly and ell=K; phi_deg is ignored.

opt.k = 0;                             % only used when task='point'
opt.ell = 25.5513726457298;             % only used when task='point'
opt.refine = true;                      % optimize K inside numerically found brackets
opt.logK_tolerance = 1e-6;              % absolute fminbnd tolerance on log K
opt.max_refine_evaluations = 120;
opt.growth_tol = 1e-6;                  % Re(sigma)>this: detected positive growth
opt.qz_tol = 1e-9;                      % denominator test for infinite constraint eigenvalues

%% 3. Checks and output options -- EDIT HERE
opt.check_N = [];                       % [] off; e.g. 128 (spectral) or 160 (staggered)
opt.check_peak = true;                  % reoptimize K as well as checking fixed-K growth
opt.check_no_growth = true;             % confirm nearly neutral Cheb results using primitives
opt.confirm_N = 160;                    % cells for that primitive check
opt.confirm_trigger = 1e-5;
opt.return_mode = true;                 % leading eigenfunction and energy diagnostics
opt.return_matrices = false;            % true saves final A/B; larger MAT file
run_reference_checks = true;            % independent fixed-wavenumber baseline checks
make_quick_plots = false;               % optional diagnostic figures, not publication Figure 3
output_directory = fullfile(here,'SI_solver_output');

%% 4. Input checks (normally do not edit below this line)
validateattributes(Ri_values,{'numeric'},{'real','finite','vector','positive'});
validateattributes(Ev_values,{'numeric'},{'real','finite','vector','nonnegative'});
Ri_values = unique(Ri_values(:).','sorted');
Ev_values = unique(Ev_values(:),'sorted');
if ~exist(output_directory,'dir'), mkdir(output_directory); end

%% 5. Optional benchmark checks -- no fitted formula enters these solves
if run_reference_checks
    % Columns: Ri, Ev, ell in manuscript units, archived Chebyshev96 growth.
    % These constants are used ONLY to check the solver, never to find peaks.
    reference = [ ...
        0.70, 1e-4, 25.5513726457298, 0.5783219138885136; ...
        0.70, 1e-3, 14.1952417050361, 0.3939999028234078; ...
        0.20, 1e-5, 135.075709407515, 1.7944680784757128; ...
        0.80, 2e-3, 11.6071539080867, 0.0264124484427052];
    reference_result = zeros(size(reference,1),6);
    fprintf('Checking four independent Chebyshev fixed-wavenumber benchmarks...\n');
    for j=1:size(reference,1)
        pt=struct('Ri',reference(j,1),'Ev',reference(j,2),'delta',.1,'Pr',1, ...
            'buoyancy_bc','neumann','horizontal_ratio',0,'U_offset',.5);
        ot=struct('task','point','method','chebyshev','N',96,'k',0, ...
            'ell',reference(j,3),'return_mode',false,'check_no_growth',false);
        test=solve_SI_instability(pt,ot);
        err=100*(test.peak.growth/reference(j,4)-1);
        reference_result(j,:)=[reference(j,:),test.peak.growth,err];
        fprintf('  Ri=%.2f  Ev=%.1e: growth=%.10f; relative difference=%+.3g%%\n', ...
            pt.Ri,pt.Ev,test.peak.growth,err);
        if abs(err)>0.05
            warning('SI:Benchmark',['Baseline growth differs by more than 0.05%%. ', ...
                'Inspect numerical conditioning and resolution before using this setup.']);
        end
    end
    Tref=array2table(reference_result,'VariableNames', ...
        {'Ri','Ev','ell','archived_growth','computed_growth','relative_difference_pct'});
    writetable(Tref,fullfile(output_directory,'fixed_wavenumber_benchmarks.csv'));
end

%% 6. Solve each chosen (Ri,Ev) case
nR=numel(Ri_values); nE=numel(Ev_values);
solutions=cell(nE,nR);
% Matrices are Ev-by-Ri, matching the earlier plotting packages.
growth_max=nan(nE,nR); sigma_complex=complex(nan(nE,nR));
sigma_star=nan(nE,nR); ell_star=nan(nE,nR); k_star=nan(nE,nR);
K_star=nan(nE,nR); lambda_over_L=nan(nE,nR);
unstable=false(nE,nR); scan_boundary=false(nE,nR);
selected_wavelength_defined=false(nE,nR);
rows=cell(nE*nR,18); q=0; start_all=tic;
for ir=1:nR
    p.Ri=Ri_values(ir);
    for ie=1:nE
        p.Ev=Ev_values(ie); start_case=tic;
        R=solve_SI_instability(p,opt); solutions{ie,ir}=R; peak=R.peak;
        growth_max(ie,ir)=peak.growth; sigma_complex(ie,ir)=peak.sigma;
        unstable(ie,ir)=peak.has_growth; scan_boundary(ie,ir)=peak.radial_boundary;
        selected_wavelength_defined(ie,ir)=peak.selected_wavelength_defined;
        if peak.selected_wavelength_defined
            sigma_star(ie,ir)=peak.growth;
            ell_star(ie,ir)=peak.ell; k_star(ie,ir)=peak.k; K_star(ie,ir)=peak.K;
            lambda_over_L(ie,ir)=peak.lambda_over_L;
        end
        q=q+1;
        rEigen=NaN; rPrimitive=NaN;
        if opt.return_mode
            rEigen=R.mode.eigen_residual; rPrimitive=R.mode.primitive_residual;
        end
        rows(q,:)={p.Ri,p.Ev,p.delta,p.Pr,opt.N,R.method_used, ...
            peak.ell,peak.k,peak.K,peak.phi_deg,peak.growth,imag(peak.sigma), ...
            peak.lambda_over_L,peak.has_growth,peak.radial_boundary,peak.status, ...
            rEigen,rPrimitive};
        fprintf('%3d/%3d Ri=%.2f Ev=%.2e  Re(sigma)=%.8f  lambda/L=%.6g  [%.1f s]\n', ...
            q,nR*nE,p.Ri,p.Ev,peak.growth,peak.lambda_over_L,toc(start_case));
        fprintf('       %s\n',peak.status);
        if ~isempty(R.notes), fprintf('       %s\n',R.notes{:}); end
        % Save after every completed case so a long sweep can be inspected.
        save(fullfile(output_directory,'SI_results.mat'),'Ri_values','Ev_values','p','opt', ...
            'solutions','growth_max','sigma_complex','sigma_star','ell_star','k_star', ...
            'K_star','lambda_over_L','unstable','scan_boundary','selected_wavelength_defined','-v7');
    end
end
summary=cell2table(rows,'VariableNames',{'Ri','Ev','delta','Pr','N','solver', ...
    'candidate_ell','candidate_k','candidate_K','candidate_phi_deg','growth_max', ...
    'sigma_imag','lambda_over_L','unstable','scan_boundary','status', ...
    'eigen_residual','primitive_residual'});
writetable(summary,fullfile(output_directory,'SI_peak_summary.csv'));

% Also save a compact plotting-compatible table for strict-SI runs.
if strcmpi(opt.task,'si')
    [Ri_matrix,Ev_matrix]=meshgrid(Ri_values,Ev_values);
    Figure3_data=table(Ri_matrix(:),Ev_matrix(:), ...
        repmat(p.delta,nE*nR,1),repmat(p.Pr,nE*nR,1),repmat(opt.N,nE*nR,1), ...
        ell_star(:),sigma_star(:),lambda_over_L(:),unstable(:),scan_boundary(:), ...
        'VariableNames',{'Ri','Ev','delta','Pr','N','ell_star','sigma_star', ...
        'lambda_over_L','unstable','scan_boundary'});
    writetable(Figure3_data,fullfile(output_directory,'SI_numerical_data.csv'));
end

%% 7. Optional quick diagnostic plots (no modification of the computed data)
if make_quick_plots && strcmpi(opt.task,'si')
    % Growth and raw wavelength are kept in distinct diagnostic windows.
    figure('Color','w'); hold on;
    for ir=1:nR
        use=Ev_values>0 & selected_wavelength_defined(:,ir);
        semilogx(Ev_values(use),sigma_star(use,ir),'-o','MarkerFaceColor','none', ...
            'DisplayName',sprintf('Ri=%.2f',Ri_values(ir)));
    end
    set(gca,'XScale','log'); xlabel('E_v'); ylabel('Nondimensional growth rate Re(\sigma)');
    legend('show','Location','best'); box on;
    figure('Color','w'); hold on;
    for ir=1:nR
        use=Ev_values>0 & selected_wavelength_defined(:,ir);
        loglog(Ev_values(use),lambda_over_L(use,ir),'o','MarkerFaceColor','none', ...
            'DisplayName',sprintf('Ri=%.2f',Ri_values(ir)));
    end
    set(gca,'XScale','log','YScale','log'); xlabel('E_v'); ylabel('\lambda_*/L');
    legend('show','Location','best'); box on;
end
fprintf('\nDone. Results saved in:\n%s\nElapsed: %.1f s\n',output_directory,toc(start_all));
fprintf('All growth rates are nondimensional. Dimensional rate = f*Re(sigma), not Re(sigma)/f.\n');
