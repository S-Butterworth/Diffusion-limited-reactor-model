function [W,Y,C] = Reactor(C,F,data)

function k = kinetics(Y,C)
    Y = Y(:);
    T = Y(end);
    k = C.A .* exp(-C.EA./((T+273.15)*C.R));
end

function Act = activity(C)
    if C.n == 1
        Act = C.a0 * exp(-C.kd * C.t);
    else
        Act = (C.a0^(1-C.n) + C.kd*(C.n-1)*C.t)^(1/(1-C.n));
    end
end

function rate = reactions(Y,C,k)
    Y = Y(:);
    Fx = Y(1:end-2);
    P = Y(end-1);

    Pi = P * Fx / sum(Fx);

    logPi = zeros(size(Pi));
    nz = Pi > 0;
    logPi(nz) = log(Pi(nz));

    contrib = C.coeff' .* logPi;
    missing = (C.coeff' ~= 0) & (Pi == 0);
    contrib(missing) = -Inf;  
    rate = C.a*k(:)' .* (exp(sum(contrib,1)));
end


function dTdW = thermodynamic(Y,C,F,rate)
    if F.thermo == false
        dTdW=0;
        return
    end

    Y = Y(:); 
    Fx = Y(1:end-2)';
    T = Y(end);

    a = C.NASA;

    Hvec = [1;T/2;(T^2)/3;(T^3)/4;(T^4)/5;1/T];
    Tvec = [1; T; T^2; T^3; T^4];

    dH = (Hvec' * a(1:6,:)) * C.R * T / 1000;
    Cp = (Tvec' * a(1:5,:))*C.R;

    deltaHr = sum(dH.*C.stoich,2);
    CpAv = sum(Fx .* Cp)/sum(Fx);
    dTdW = -1000 * sum(deltaHr .* rate') / (sum(Fx) * CpAv);
end

function dPdW = pressure(Y,C,F)
    if F.pressure == false
        dPdW=0;
        return
    end

    Y = Y(:);
    Fx = Y(1:end-2);
    P = Y(end-1);
    T = Y(end);

    n2 = - C.G * (150 * (1 - C.vf) * C.mu/C.dp + 1.75 * C.G) * C.PT * sum(Fx) * T;
    d2 = C.dens * C.dp * C.vf^3 * C.pc * C.CSA * P * C.F0 * C.T0;

    dPdW = n2 / d2;
end

function eta = lookupEta(Y,C,F,data)
    if F.diffusion == false
        eta=1;
        return
    end
    eta = interp2(data.t_grid, data.T_grid, data.eta_table, C.t, Y(end), 'linear');
    if isnan(eta)
        error('Requested (T=%.1f, t=%.1f) falls outside the validated table range.',C.t, Y(end));
    end
end

function dydt = pfr(~,Y,C,F,data)
    Y = Y(:); 

    k = kinetics(Y,C);

    eta = lookupEta(Y,C,F,data);

    rate = eta * reactions(Y,C,k);
    dFdW = sum(C.stoich' .* rate,2);

    dTdW = thermodynamic(Y,C,F,rate);

    dPdW = pressure(Y,C,F);

    dydt = [dFdW;dPdW;dTdW];
end

if F.activity == false
    C.a=1;
else
    C.a = activity(C);
end

Y0 = [C.F,C.PT,C.T0];
[W,Y] = ode45(@(W,Y) pfr(W,Y,C,F,data),[0,C.cat],Y0);

end