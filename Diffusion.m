function [T0_steps, eta_history, rc_history,sol] = Diffusion(T0_steps, C)
    C.a = activity(C);
    S = [0, 0, 0, 0;
         0,-2, 0, 0;
         0, 0, 0, 0;
         0, 0, 0,-2];
    options1 = bvpset('SingularTerm', S, 'RelTol', C.RelTol1);
    options2 = bvpset('RelTol', C.RelTol2);

    usingFreeBoundary = false;

    mesh = [linspace(0,1,1000)];
    solinit = bvpinit(mesh,[0.5;0;0.5;0]);

    n = length(T0_steps);
    eta_history = nan(1,n);
    rc_history  = nan(1,n);
    delta_reg = 1e-6;

    for i = 1:n
        C.T = T0_steps(i) + 273.15;
        if ~usingFreeBoundary
            sol = bvp4c(@(r,y) diffODE(r,y,C), @(ya,yb) diffBC(ya,yb), solinit, options1);
            solinit = sol;

            if min(sol.y(1,:)) < C.threshold
                idx = find(sol.y(1,:) < C.threshold, 1, 'last');
                rc_guess = sol.x(idx);
                if rc_guess >= C.minRcToSwitch
                    trial.x = sol.x;
                    trial.y = sol.y;
                    trial.parameters = rc_guess;

                    testsol = bvp4c(@(s,y,rc) fbODE(s,y,rc,C), @(ya,yb,rc) fbBC(ya,yb,rc,C), trial, options2);
                    if testsol.parameters > 0
                        sol = testsol;
                        usingFreeBoundary = true;
                        solinit = sol;
                    end
                end
            end
        else
            sol = bvp4c(@(s,y,rc) fbODE(s,y,rc,C), @(ya,yb,rc) fbBC(ya,yb,rc,C), solinit, options2);
            solinit = sol;
        end

        if usingFreeBoundary
            rc = sol.parameters;
            r_active = rc + (1-rc)*sol.x;
            CA_floor = 0.5*(sol.y(1,:) + sqrt(sol.y(1,:).^2 + delta_reg^2));
            CB_floor = 0.5*(sol.y(3,:) + sqrt(sol.y(3,:).^2 + delta_reg^2));
            localRate = CA_floor.^(1/3) .* CB_floor.^(2/3);
            numerator = trapz(sol.x, localRate .* r_active.^2) * (1-rc);
            rc_history(i) = rc;
        else
            CA_floor = 0.5*(sol.y(1,:) + sqrt(sol.y(1,:).^2 + delta_reg^2));
            CB_floor = 0.5*(sol.y(3,:) + sqrt(sol.y(3,:).^2 + delta_reg^2));
            localRate = CA_floor.^(1/3) .* CB_floor.^(2/3);
            numerator = trapz(sol.x, localRate .* sol.x.^2);
        end
        denominator = 1/3;
        eta_history(i) = numerator/denominator;
    end
end

function k = kinetics(T,C)
    k = C.A .* exp(-C.EA./(T*C.R));
end

function Act = activity(C)
    if C.n == 1
        Act = C.a0 * exp(-C.kd * C.t);
    else
        Act = (C.a0^(1-C.n) + C.kd*(C.n-1)*C.t)^(1/(1-C.n));
    end
end

function dydx = fbODE(s,y,rc,C)
    PA = C.x(1) * C.PT;
    PB = C.x(2) * C.PT;
    CAs = PA / (C.R*C.T);
    CBs = PB / (C.R*C.T);

    k = kinetics(C.T,C);

    phiA2 = C.a * ((C.dp/2)^2 * k * C.pc * C.R * C.T * (CAs^(1/3 -1)) * (CBs^(2/3))) / C.D(1);
    phiB2 = C.a * 0.5 * ((C.dp/2)^2 * k * C.pc * C.R * C.T * (CAs^(1/3)) * (CBs^(2/3 -1))) / C.D(2);

    r = rc + (1-rc)*s;
    delta_reg  = 1e-6;
    CA_floor = 0.5*(y(1) + sqrt(y(1)^2 + delta_reg ^2));
    CB_floor = 0.5*(y(3) + sqrt(y(3)^2 + delta_reg ^2));
    rate = CA_floor^(1/3) * CB_floor^(2/3);

    dydx = [(1-rc)*y(2);
            (1-rc)*(phiA2*rate - (2/r)*y(2));
            (1-rc)*y(4);
            (1-rc)*(phiB2*rate - (2/r)*y(4))];
end

function dydx = diffODE(~,y,C)
    PA = C.x(1) * C.PT;
    PB = C.x(2) * C.PT;
    CAs = PA / (C.R*C.T);
    CBs = PB / (C.R*C.T);

    k = kinetics(C.T,C);

    phiA2 = C.a * ((C.dp/2)^2 * k * C.pc * C.R * C.T * (CAs^(1/3 -1)) * (CBs^(2/3))) / C.D(1);
    phiB2 = C.a * 0.5 * ((C.dp/2)^2 * k * C.pc * C.R * C.T * (CAs^(1/3)) * (CBs^(2/3 -1))) / C.D(2);

    delta_reg = 1e-6;
    CA_floor = 0.5*(y(1) + sqrt(y(1)^2 + delta_reg^2));
    CB_floor = 0.5*(y(3) + sqrt(y(3)^2 + delta_reg^2));
    rate = CA_floor^(1/3) * CB_floor^(2/3);
    
    dydx = [y(2); phiA2*rate;
            y(4); phiB2*rate];
end

function res = fbBC(ya, yb, ~, ~) 
res = [ya(1); ya(2); ya(4);   
        yb(1)-1; yb(3)-1];      
end

function res = diffBC(ya,yb)
    res = [ya(2); ya(4); 
        yb(1)-1; yb(3)-1];
end