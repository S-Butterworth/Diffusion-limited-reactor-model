%% Importing inital conditions, constants and diffusion data files
[C,F,data] = Initial_Conditions_and_Constants();

%% Runs reactor model Converging at X
function err = catErrorFn(cat_trial, C, F, data)
    Clocal = C;
    Clocal.cat = cat_trial;
    [~, Y, ~] = Reactor(Clocal, F, data);
    achieved_X = (Y(1,1) - Y(end,1)) / Y(1,1);
    err = achieved_X - Clocal.X;
end

function [W,Y,C,converged] = converge(C,F,data)
    cat_max = 1e6;

    warnState = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
    cleanupObj = onCleanup(@() warning(warnState));

    Clocal = C; Clocal.cat = cat_max;
    [~, Y_max, ~] = Reactor(Clocal, F, data);
    achievable_X = (Y_max(1,1) - Y_max(end,1)) / Y_max(1,1);

    if achievable_X < C.X
        fprintf('Target X=%.2f unreachable even at cat=%.2e; plateaus at %.4f\n', ...
                 C.X, cat_max, achievable_X);
        C.cat = cat_max;
        [W,Y,C] = Reactor(C,F,data);
        converged = false;
        return
    end

    opts = optimset('TolX', 1e-8);
    cat_solution = fzero(@(cat) catErrorFn(cat, C, F, data), [C.cat*0.01, cat_max], opts);
    C.cat = cat_solution;
    [W,Y,C] = Reactor(C,F,data);
    converged = true;
end

[W,Y,C] = converge(C,F,data);

%% Displays Reactor set up settings and results
if F.stats == true
    fprintf('=== Reactor Run Settings ===\n');
    fprintf('Inlet temperature (T0):      %.1f C\n', C.T0);
    fprintf('Inlet pressure (PT):         %.3e Pa\n', C.PT);
    fprintf('Catalyst time-in-use (t):    %.1f days\n', C.t);
    fprintf('Target/achieved conversion:  %.2f\n', C.X);
    fprintf('Reactor length (L):          %.2f m\n', W(end) /(C.CSA * C.pc * (1-C.vf)));
    fprintf('Catalyst weight (kg):        %.2f kg\n', C.NoReactors*W(end));
    fprintf('Pressure drop:               %.3e Pa\n', Y(1,5)-Y(end,5));
    fprintf('Feed mole fractions (xA,xB): %.4f, %.4f\n', C.x(1), C.x(2));

    fprintf('\n=== Kinetics & Deactivation ===\n');
    fprintf('Rate constant (k, ref T):    %.4e\n', C.k);
    fprintf('Activation energy (EA):      %.1f J/mol\n', C.EA);
    fprintf('Deactivation constant (kd):  %.4e\n', C.kd);
    fprintf('Deactivation order (n):      %.1f\n', C.n);
    fprintf('Activity this run (Act):     %.4f\n', C.a);

    fprintf('\n=== Catalyst Pellet ===\n');
    fprintf('Pellet diameter (dp):        %.4e m\n', C.dp);
    fprintf('Catalyst density (pc):       %.1f kg/m3\n', C.pc);
    fprintf('Diffusivities (DA, DB):      %.4e, %.4e m2/s\n', C.D(1), C.D(2));
end

%% File directory for figures
iptDir = fileparts(mfilename('fullpath'));
outDir = fullfile(iptDir, 'figures');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% Plotting mole fraction vs catalyst weight
mf = Y(:,1:3)./sum(Y(:,1:4),2);

f = figure('Visible', 'on');
theme(f, "light");
ax = axes(f);

plot(ax, W /(C.CSA * C.pc * (1-C.vf)), mf, 'LineWidth', 1.2)
grid(ax, 'on')
xlim(ax,[0,W(end) /(C.CSA * C.pc * (1-C.vf))])
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.05;

legend(ax, ["A", "B", "C"], 'Location', 'northwest', 'Box', 'off')
title(ax, "Mole Fraction vs Reactor Length")
xlabel(ax, "Reactor Length (m)")
ylabel(ax, "Mole Fraction")

exportgraphics(f, fullfile(outDir, 'mole_fraction.png'), ...
    'Resolution', 600, 'BackgroundColor', 'white', 'Padding', 60);

%% Plotting Pressure and temperature against catalyst weight
f = figure('Visible', 'on');
theme(f, "light");
ax = axes(f);

plot(ax, W /(C.CSA * C.pc * (1-C.vf)),Y(:,5)/1E5, 'LineWidth', 1.2)
ylabel(ax, "Pressure (bar)")
yyaxis right
plot(ax, W /(C.CSA * C.pc * (1-C.vf)),Y(:,6), 'LineWidth', 1.2)
grid(ax, 'on')
xlim(ax,[0,W(end) /(C.CSA * C.pc * (1-C.vf))])
ax.YAxis(2).Color = "k"; 
ax.YAxis(2).TickLabelFormat = ' %g';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.05;

legend(ax, ["Pressure", "Temperature"], 'Location', 'southeast', 'Box', 'off')
title(ax, "Temperature and Pressure vs Reactor Length")
xlabel(ax, "Reactor Length (m)")
ylabel(ax, "Temperature (°C)")

exportgraphics(f, fullfile(outDir, 'temperature_and_pressure.png'), ...
     'Resolution', 600, 'BackgroundColor', 'white', 'Padding', 60);

%% Plotting mole fraction vs reactor length
mf = Y(:,1:3)./sum(Y(:,1:4),2);

f = figure('Visible', 'on');
theme(f, "light");
ax = axes(f);

plot(ax, W /(C.CSA * C.pc * (1-C.vf)), (Y(1,1)-Y(:,1))/Y(1,1), 'LineWidth', 1.2)
grid(ax, 'on')
xlim(ax,[0,W(end) /(C.CSA * C.pc * (1-C.vf))])
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.05;

title(ax, "Conversion vs Reactor Length")
xlabel(ax, "Reactor Length (m)")
ylabel(ax, "Conversion")

exportgraphics(f, fullfile(outDir, 'conversion.png'), ...
    'Resolution', 600, 'BackgroundColor', 'white', 'Padding', 60);

%% Conversion as catalyst decays
function Convals = decay_curve(C,F,data)
    Convals = zeros(1,101);
    for t = 0:100
        C.t = t;
        [~,Y,C] = Reactor(C,F,data);
        Convals(t+1) = ((Y(1,1)-Y(end,1))/Y(1,1));
    end 
end
Convals = decay_curve(C,F,data);

Flocal = F;
Flocal.diffusion = false;
[Wlocal,~,Clocal] = converge(C,Flocal,data);
ConvalsNoDiffusion = decay_curve(Clocal,Flocal,data);

f = figure('Visible', 'on');
theme(f, "light");
ax = axes(f);

t = linspace(0,100,101);
plot(ax, t, Convals, 'LineWidth', 1.2)
hold on;
plot(ax, t, ConvalsNoDiffusion, 'LineWidth', 1.2)
grid(ax, 'on')
xlim(ax,[0,100])
legend(ax, ["Diffusion:ON", "Diffusion:OFF"], 'Location', 'northeast', 'Box', 'off')
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.05;

title(ax, "Conversion vs time on stream")
xlabel(ax, "Time on stream (days)")
ylabel(ax, "Conversion")
hold off;

exportgraphics(f, fullfile(outDir, 'conversion_vs_TOS.png'), ...
    'Resolution', 600, 'BackgroundColor', 'white', 'Padding', 60);

%% Required Inlet temperature

error = zeros(1,101);
increased_temp = zeros(1,101);
T0 = C.T0;
function err = convErrorFn(T0_trial, C, F, data)
    Clocal = C;
    Clocal.T0 = T0_trial;
    [~, Y, ~] = Reactor(Clocal, F, data);
    achieved_X = (Y(1,1) - Y(end,1)) / Y(1,1);
    err = achieved_X - Clocal.X;
end

for t = 0:100
    C.t = t;

    T0_solution = fzero(@(T0) convErrorFn(T0, C, F, data), C.T0);
    C.T0 = T0_solution;
    [W, Y, C] = Reactor(C, F, data);

    increased_temp(t+1) = C.T0 - T0;
    achieved_X = (Y(1,1) - Y(end,1)) / Y(1,1);
    error(t+1) = 100*abs(C.X - achieved_X);
end

f = figure('Visible', 'on');
theme(f, "light");
ax = axes(f);

t = linspace(0,100,101);
plot(ax, t, increased_temp, 'LineWidth', 1.2)
ylabel(ax, "Feed temperature increase (°C)")
grid(ax, 'on')
xlim(ax,[0,100])
ylim(ax,[0,inf])
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.05;

title(ax, "Temperature increase for stable conversion rates")
xlabel(ax, "Time on stream (days)")

exportgraphics(f, fullfile(outDir, 'Temperature_increase_for_constant_conversion.png'), ...
    'Resolution', 600, 'BackgroundColor', 'white', 'Padding', 60);