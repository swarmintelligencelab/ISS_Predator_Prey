%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project : Input-to-State Stability of Tracking Dynamics in a Predator-Prey Feedback Control Model
%Author  : Deze Liu, Daniel Burbano (db1359@soe.rutgers.edu)
%Lab     : The Swarm Intelligence Lab
%Date    : 09/25/2026
%Description :
% This MATLAB code simulates the predator–prey dynamics and evaluates the theoretical tracking-error bounds. 
% It compares the simulated bearing error with the bounds derived in Theorem 1 and Corollary 1, 
% and performs parameter sweeps to investigate how predator and prey control parameters affect 
% the realized tracking error and its theoretical bound.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear;
close all;


%% Initial setting
cfg.T       = 30;
cfg.maxStep = 1/180;
cfg.relTol  = 1e-9;
cfg.absTol  = 1e-11;

cfg.r_star = 2.0;

cfg.gammaGate = 100;

cfg.kappa_d_ref = 42.64;
cfg.mu_ref      = 0.73;
cfg.tau_d       = 1/sqrt(cfg.kappa_d_ref*cfg.mu_ref);

cfg.epsilon_d = 1.0;

cfg.nSweep = 21;
cfg.sweepLowerFactor = 0.5;
cfg.sweepUpperFactor = 1.5;

cfg.tol = 1e-6;

% Figure output
cfg.saveResolution = 300;
cfg.figureSizeIn = [7.5 5.5];
cfg.lineWidth = 2.6;
cfg.axisLineWidth = 2.0;
cfg.startMarkerSize = 16;
cfg.fontSize = 28;
cfg.labelFontSize = 32;
cfg.axesPosition = [0.21 0.23 0.75 0.64];
cfg.yPaddingFraction = 0.06;
cfg.xPaddingFraction = 0.02;
cfg.nSweepTicks = 3;
scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    outputDir = pwd;
else
    outputDir = fileparts(scriptFullPath);
end

cfg.z0 = [sqrt(30); 0; 0.01; 0; ...   % fish:    x,y,theta,omega
          0;        0; 0.01; 0];      % dolphin: x,y,theta,omega


%% Parameters
before.name    = 'Before';
before.alpha_f = 23.95;
before.kappa_f = 56.33;
before.delta   = 78.33;
before.d_s     = 4.77;
before.gamma   = cfg.gammaGate;

before.v_f     = 8.95;
before.v_d     = 10.01;

before.alpha_d = 21.27;
before.kappa_d = 42.64;
before.mu      = 0.73;

R = runCase(before,cfg);

fprintf('\nBefore Learning\n');
fprintf('  initial r(0)            = %.4f m\n',R.r(1));
fprintf('  min r(t)                = %.4f m\n',min(R.r));
fprintf('  max |e_d(t)|            = %.4f rad\n',max(abs(R.e_d)));
fprintf('  max Theorem 1 bound     = %.4f rad\n',max(R.Bconv));
fprintf('  max Corollary 1 bound   = %.4f rad\n',max(R.Bparam));
fprintf('  max |zeta_d(t)|         = %.4f rad/s^2\n',max(abs(R.zeta_d)));
fprintf('  Corollary zeta bound    = %.4f rad/s^2\n',R.zetaBar);
fprintf('  e_d remains in (-pi,pi) = %d\n',R.angleValid);
fprintf('  r(t) >= r_star          = %d\n',min(R.r) >= cfg.r_star-cfg.tol);


%% Parameter sweeps
kappaDVals = linspace(cfg.sweepLowerFactor*before.kappa_d,...
                      cfg.sweepUpperFactor*before.kappa_d,...
                      cfg.nSweep);
alphaDVals = linspace(cfg.sweepLowerFactor*before.alpha_d,...
                      cfg.sweepUpperFactor*before.alpha_d,...
                      cfg.nSweep);
muVals     = linspace(cfg.sweepLowerFactor*before.mu,...
                      cfg.sweepUpperFactor*before.mu,cfg.nSweep);
kappaFVals = linspace(cfg.sweepLowerFactor*before.kappa_f,...
                      cfg.sweepUpperFactor*before.kappa_f,...
                      cfg.nSweep);

SkappaD = sweepParameter(before,'kappa_d',kappaDVals,cfg);
SalphaD = sweepParameter(before,'alpha_d',alphaDVals,cfg);
Smu     = sweepParameter(before,'mu',muVals,cfg);
SkappaF = sweepParameter(before,'kappa_f',kappaFVals,cfg);


%% Plot figures
% (a) Before-Learning planar trajectories.
fig = newPaperFigure('Figure 2a - Before trajectories',cfg);
ax = axes(fig);
fishPath = plot(ax,R.z(:,1),R.z(:,2),'LineWidth',cfg.lineWidth);
hold(ax,'on');
dolphinPath = plot(ax,R.z(:,5),R.z(:,6),'LineWidth',cfg.lineWidth);
plot(ax,R.z(1,1),R.z(1,2),'*',...
     'Color',[0 0.4470 0.7410],...
     'MarkerSize',cfg.startMarkerSize,...
     'LineWidth',cfg.lineWidth,...
     'HandleVisibility','off');
plot(ax,R.z(1,5),R.z(1,6),'*',...
     'Color',[1 0 0],...
     'MarkerSize',cfg.startMarkerSize,...
     'LineWidth',cfg.lineWidth,...
     'HandleVisibility','off');
axis(ax,'equal');
xlim(ax,[-160 100]);
ylim(ax,[-150 15]);
grid(ax,'on');
xlabel(ax,'x [m]');
ylabel(ax,'y [m]');
legend(ax,[fishPath dolphinPath],...
       '$\mathrm{Prey}$','$\mathrm{Predator}$',...
       'Interpreter','latex','Location','best');
applyPaperAxesStyle(ax,cfg);
savePaperFigure(fig,outputDir,'Fig2a.png',...
                cfg.saveResolution);

% (b) Before-Learning tracking error and Theorem 1 bound.
fig = newPaperFigure('Figure 2b - Before tracking bound',cfg);
ax = axes(fig);
plot(ax,R.t,abs(R.e_d),'LineWidth',cfg.lineWidth); hold(ax,'on');
plot(ax,R.t,R.Bconv,'--','LineWidth',cfg.lineWidth);
grid(ax,'on');
xlabel(ax,'Time [s]');
ylabel(ax,'Bearing error [rad]');
applyPaddedYLimits(ax,[abs(R.e_d(:)); R.Bconv(:)],...
                   cfg.yPaddingFraction,true);
legend(ax,'$|e_d(t)|$','$B_{\mathrm{conv}}(t)$',...
       'Interpreter','latex','Location','best');
applyPaperAxesStyle(ax,cfg);
savePaperFigure(fig,outputDir,'Fig2b.png',...
                cfg.saveResolution);

% (c) Before-Learning LOS forcing and Corollary bound.
fig = newPaperFigure('Figure 2c - Before LOS forcing',cfg);
ax = axes(fig);
plot(ax,R.t,abs(R.zeta_d),'LineWidth',cfg.lineWidth); hold(ax,'on');
yline(ax,R.zetaBar,'--','LineWidth',cfg.lineWidth);
grid(ax,'on');
xlabel(ax,'Time [s]');
ylabel(ax,'|\zeta_d| [rad/s^2]');
applyPaddedYLimits(ax,[abs(R.zeta_d(:)); R.zetaBar],...
                   cfg.yPaddingFraction,true);
legend(ax,'$|\zeta_d(t)|$','$\bar{\zeta}_d$',...
       'Interpreter','latex','Location','best');
applyPaperAxesStyle(ax,cfg);
savePaperFigure(fig,outputDir,'Fig2c.png',...
                cfg.saveResolution);

% (d) Before-Learning tracking error and Corollary 1 parameter-only bound.
fig = newPaperFigure('Figure 2d - Corollary 1 tracking bound',cfg);
ax = axes(fig);
plot(ax,R.t,abs(R.e_d),'LineWidth',cfg.lineWidth); hold(ax,'on');
plot(ax,R.t,R.Bparam,'--','LineWidth',cfg.lineWidth);
grid(ax,'on');
xlabel(ax,'Time [s]');
ylabel(ax,'Bearing error [rad]');
applyPaddedYLimits(ax,[abs(R.e_d(:)); R.Bparam(:)],...
                   cfg.yPaddingFraction,true);
legend(ax,'$|e_d(t)|$','$B_{\mathrm{param}}(t)$',...
       'Interpreter','latex','Location','east');
applyPaperAxesStyle(ax,cfg);
savePaperFigure(fig,outputDir,'Fig2d.png',...
                cfg.saveResolution);

% (e)-(h) Parameter sweeps.
fig = newPaperFigure('Figure 2e - kappa_d sweep',cfg);
plotSweep(axes(fig),SkappaD,'\kappa_d',cfg);
savePaperFigure(fig,outputDir,'Fig2e.png',...
                cfg.saveResolution);

fig = newPaperFigure('Figure 2f - alpha_d sweep',cfg);
plotSweep(axes(fig),SalphaD,'\alpha_d [s^{-1}]',cfg);
savePaperFigure(fig,outputDir,'Fig2f.png',...
                cfg.saveResolution);

fig = newPaperFigure('Figure 2g - mu sweep',cfg);
plotSweep(axes(fig),Smu,'\mu [rad^{-1}]',cfg);
savePaperFigure(fig,outputDir,'Fig2g.png',...
                cfg.saveResolution);

fig = newPaperFigure('Figure 2h - kappa_f sweep',cfg);
plotSweep(axes(fig),SkappaF,'\kappa_f',cfg);
savePaperFigure(fig,outputDir,'Fig2h.png',...
                cfg.saveResolution);

fprintf('\nSaved eight figures to:\n  %s\n',outputDir);





function R = runCase(p,cfg)
    opts = odeset('RelTol',cfg.relTol,'AbsTol',cfg.absTol,...
                  'MaxStep',cfg.maxStep);

    [t,z] = ode45(@(t,z) modelODE(t,z,p),[0 cfg.T],cfg.z0,opts);

    n = numel(t);

    r      = zeros(n,1);
    e_f    = zeros(n,1);
    e_d    = zeros(n,1);
    q      = zeros(n,1);
    qdot   = zeros(n,1);
    edot   = zeros(n,1);
    zeta_d = zeros(n,1);

    for k = 1:n
        [r(k),e_f(k),e_d(k)] = geometry(z(k,:).');

        omega_f = z(k,4);
        omega_d = z(k,8);

        q(k) = (p.v_d*sin(e_d(k)) + p.v_f*sin(e_f(k)))/r(k);
        edot(k) = q(k)-omega_d;

        G = p.v_d*cos(e_d(k)) + p.v_f*cos(e_f(k));
        J = p.v_d*cos(e_d(k))*omega_d + ...
            p.v_f*cos(e_f(k))*omega_f;

        qdot(k) = (2*q(k)*G-J)/r(k);
        zeta_d(k) = qdot(k)+p.alpha_d*q(k);
    end

    % Theorem 1 constants.
    eps_d = cfg.epsilon_d;
    tau_d = cfg.tau_d;

    if eps_d <= 0 || eps_d >= p.alpha_d
        error('epsilon_d must satisfy 0 < epsilon_d < alpha_d.');
    end

    ell_d = tanh(p.mu*pi)/pi;

    aPlus = tau_d^2*(p.kappa_d*p.mu + eps_d*p.alpha_d);
    b_d   = eps_d*tau_d;

    lambdaOver = 0.5*(...
        1 + aPlus + sqrt((aPlus-1)^2 + 4*b_d^2));

    nu_d = min(p.alpha_d-eps_d,...
               eps_d*p.kappa_d*ell_d*tau_d^2);

    lambda_d = nu_d/lambdaOver;

    c_d = tau_d^2/(2*(p.alpha_d-eps_d)) + ...
          eps_d*tau_d^2/(2*p.kappa_d*ell_d);

    beta_d = tau_d^2*(...
        p.kappa_d*ell_d + eps_d*(p.alpha_d-eps_d));

    % Scaled ISS-Lyapunov function W_d in completed-square form.
    Wd = tau_d^2*(...
        0.5*(edot + eps_d*e_d).^2 + ...
        0.5*eps_d*(p.alpha_d-eps_d)*e_d.^2 + ...
        (p.kappa_d/p.mu)*log(cosh(p.mu*e_d)));

    elapsed = t-t(1);
    weightedInput = exp(lambda_d*elapsed).*zeta_d.^2;
    memory = exp(-lambda_d*elapsed).*cumtrapz(t,weightedInput);

    Wbound = exp(-lambda_d*elapsed)*Wd(1) + c_d*memory;
    Bconv = sqrt(max(2*Wbound/beta_d,0));

    % Corollary parameter-only forcing bound.
    Vsum = p.v_d+p.v_f;

    omegaBar_d = max(abs(cfg.z0(8)),p.kappa_d/p.alpha_d);
    omegaBar_f = max(abs(cfg.z0(4)),p.kappa_f*pi/p.alpha_f);

    zetaBar = Vsum^2/cfg.r_star^2 + ...
        (p.v_d*omegaBar_d + p.v_f*omegaBar_f)/cfg.r_star + ...
        p.alpha_d*Vsum/cfg.r_star;

    % Corollary 1 parameter-only tracking-error bound.
    WparamBound = exp(-lambda_d*elapsed).*Wd(1) + ...
        (c_d/lambda_d)*(1-exp(-lambda_d*elapsed))*zetaBar^2;
    Bparam = sqrt(max(2*WparamBound/beta_d,0));

    angleValid = max(abs(e_d)) < pi-cfg.tol && ...
                 ~any(abs(diff(e_d)) > pi);

    R.name       = p.name;
    R.p          = p;
    R.t          = t;
    R.z          = z;
    R.r          = r;
    R.e_f        = e_f;
    R.e_d        = e_d;
    R.q          = q;
    R.qdot       = qdot;
    R.edot       = edot;
    R.zeta_d     = zeta_d;
    R.lambda_d   = lambda_d;
    R.c_d        = c_d;
    R.beta_d     = beta_d;
    R.Wd         = Wd;
    R.Bconv      = Bconv;
    R.Bparam     = Bparam;
    R.zetaBar    = zetaBar;
    R.angleValid = angleValid;
end


function dz = modelODE(~,z,p)
    theta_f = z(3);
    omega_f = z(4);

    theta_d = z(7);
    omega_d = z(8);

    [r,e_f,e_d] = geometry(z);

    Rgate = 0.5*(1-tanh(0.5*p.delta*(r-p.d_s)));
    Hgate = 0.5*(1-tanh(p.gamma*cos(e_f)));
    c = Rgate*Hgate;

    u_f = p.kappa_f*c*e_f;
    u_d = p.kappa_d*tanh(p.mu*e_d);

    dz = zeros(8,1);

    dz(1) = p.v_f*cos(theta_f);
    dz(2) = p.v_f*sin(theta_f);
    dz(3) = omega_f;
    dz(4) = -p.alpha_f*omega_f + u_f;

    dz(5) = p.v_d*cos(theta_d);
    dz(6) = p.v_d*sin(theta_d);
    dz(7) = omega_d;
    dz(8) = -p.alpha_d*omega_d + u_d;
end


function [r,e_f,e_d] = geometry(z)
    x_f = z(1);
    y_f = z(2);
    theta_f = z(3);

    x_d = z(5);
    y_d = z(6);
    theta_d = z(7);

    dx = x_d-x_f;
    dy = y_d-y_f;

    r = hypot(dx,dy);
    psi = atan2(dy,dx);

    e_f = atan2(sin(psi-theta_f),cos(psi-theta_f));
    e_d = atan2(sin(psi+pi-theta_d),cos(psi+pi-theta_d));
end


function S = sweepParameter(base,fieldName,values,cfg)
    n = numel(values);

    actualMax  = nan(n,1);
    theoremMax = nan(n,1);
    valid      = false(n,1);

    for j = 1:n
        p = base;
        p.(fieldName) = values(j);
        p.name = 'Sweep';

        R = runCase(p,cfg);

        actualMax(j) = max(abs(R.e_d));
        valid(j) = R.angleValid;

        if valid(j)
            theoremMax(j) = max(R.Bconv);
        end
    end

    S.x          = values(:);
    S.actualMax  = actualMax;
    S.theoremMax = theoremMax;
    S.valid      = valid;
end


function plotSweep(ax,S,xText,cfg)
    plot(ax,S.x,S.actualMax,'-','LineWidth',cfg.lineWidth); hold(ax,'on');
    plot(ax,S.x,S.theoremMax,'--','LineWidth',cfg.lineWidth);

    grid(ax,'on');
    xlabel(ax,xText);
    ylabel(ax,'Peak error [rad]');

    xMin = min(S.x);
    xMax = max(S.x);
    xSpan = xMax-xMin;
    xPad = cfg.xPaddingFraction*xSpan;
    xTickValues = linspace(xMin,xMax,cfg.nSweepTicks);

    xlim(ax,[xMin-xPad,xMax+xPad]);
    xticks(ax,xTickValues);
    if xSpan < 1
        xticklabels(ax,compose('%.2f',xTickValues));
    else
        xticklabels(ax,compose('%.1f',xTickValues));
    end

    applyPaddedYLimits(ax,[S.actualMax(:); S.theoremMax(:)],...
                       cfg.yPaddingFraction,false);

    legend(ax,'$\max |e_d(t)|$',...
              '$\max B_{\mathrm{conv}}(t)$',...
              'Interpreter','latex',...
              'Location','best');
    applyPaperAxesStyle(ax,cfg);
end


function fig = newPaperFigure(figureName,cfg)
    fig = figure('Color','w',...
                 'Units','inches',...
                 'Position',[1 1 cfg.figureSizeIn],...
                 'PaperPositionMode','auto',...
                 'Name',figureName,...
                 'ToolBar','figure',...
                 'MenuBar','figure',...
                 'NumberTitle','off');

    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 cfg.figureSizeIn];
    fig.PaperSize = cfg.figureSizeIn;
end


function applyPaddedYLimits(ax,yValues,paddingFraction,includeZero)
    yValues = yValues(isfinite(yValues));
    if isempty(yValues)
        return;
    end

    yMin = min(yValues);
    yMax = max(yValues);
    ySpan = yMax-yMin;
    if ySpan <= eps(max(1,max(abs(yValues))))
        ySpan = max(1,abs(yMax));
    end

    yPad = paddingFraction*ySpan;
    lowerLimit = yMin-yPad;
    upperLimit = yMax+yPad;

    if includeZero && yMin >= 0
        lowerLimit = 0;
    end

    ylim(ax,[lowerLimit upperLimit]);
end


function applyPaperAxesStyle(ax,cfg)
    disableDefaultInteractivity(ax);
    ax.Toolbar.Visible = 'off';
    ax.Units = 'normalized';
    ax.PositionConstraint = 'innerposition';
    ax.Position = cfg.axesPosition;
    ax.FontSize = cfg.fontSize;
    ax.FontWeight = 'normal';
    ax.LineWidth = cfg.axisLineWidth;
    ax.TickDir = 'out';

    ax.XLabel.FontSize = cfg.labelFontSize;
    ax.YLabel.FontSize = cfg.labelFontSize;
    ax.XLabel.FontWeight = 'normal';
    ax.YLabel.FontWeight = 'normal';

    legends = findobj(ancestor(ax,'figure'),'Type','Legend');
    set(legends,'FontSize',cfg.fontSize);
end


function savePaperFigure(fig,outputDir,fileName,resolution)
    outputPath = fullfile(outputDir,fileName);
    print(fig,outputPath,'-dpng',sprintf('-r%d',resolution));
end
