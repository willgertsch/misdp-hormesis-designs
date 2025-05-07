% Find designs for Pozuelo-Campos, Casero-Alonso, Amo-Salas et al (2023)
% K=100
% max_dose=1
% theta = [3.086, 4.844,-13.113]
% N=30
% type = 'exact'
% obj = "D"
% time = 120
function [x,w] = find_design_hormesis_count(K, max_dose, theta, N, type, obj, p, alpha, time)
% Find designs using YALMIP

% compute design matrix
%dose = linspace(0.0, max_dose, K);
dose = 0.0:(max_dose/K):max_dose;
eta = theta(1) + theta(2) * dose + theta(3) * dose.^2;
v = exp(eta);
grad = @(x) [
    ones(1, K+1); 
    x; 
    x.^2
    ]';
F = grad(dose);

% set up optimization problem in YALMIP
% different versions for approximate and exact
if (type == "exact")

    n = intvar(K+1, 1); % optimize integer allocations
    M = eye(3) * 1e-11; % add 10^-11 to avoid singularities
    for i = 1:K+1
        M = M + n(i) * F(i, :)' * F(i, :) * v(i);
    end

    constraints = [n >= 0, sum(n) == N];
    options = sdpsettings('solver','bnb','bnb.solver','mosek', 'bnb.maxtime', time);

elseif (type == "approx")

    w = sdpvar(K+1, 1); % optimize weights
    M = eye(3) * 1e-11; % add 10^-11 to avoid singularities
    for i = 1:K+1
        M = M + w(i) * F(i, :)' * F(i, :) * v(i);
    end
    
    constraints = [w>=0, sum(w) == 1];
    options = sdpsettings('solver', 'mosek');
end

if (obj == "D")
    objective = -logdet(M);
elseif (obj == "theta0")

    c = [1;0;0];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;
elseif (obj == "theta1")

    c = [0;1;0];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;
elseif (obj == "theta2")

    c = [0;0;1];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;
elseif (obj == "zep")
    c = [0; -1/theta(3); theta(2)/theta(3)^2];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;

elseif (obj == "rauc")

    c = c_rauc(theta);
    
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;

elseif (obj == "mrd")
    c = [0; -1/(2*theta(3)); theta(2)/(2*theta(3)^2)];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;

elseif (obj == "rauc+mrd")

    c1 = c_rauc(theta);
    c2 = [0; -1/(2*theta(3)); theta(2)/(2*theta(3)^2)];
    % epigraph + Schur complement for both
    z1 = sdpvar(1);
    z2 = sdpvar(1);
    Mschur1 = [M c1 ; c1' z1];
    Mschur2 = [M c2 ; c2' z2];
    constraints = [constraints, Mschur1 >= 0, Mschur2 >= 0];
    objective = z1 + z2;
    % equivalent to compound objective with alpha = 0.5?


elseif (obj == "rauc+zep")

    c1 = c_rauc(theta);
    c2 = [0; -1/theta(3); theta(2)/theta(3)^2];
    z1 = sdpvar(1);
    z2 = sdpvar(1);
    Mschur1 = [M c1 ; c1' z1];
    Mschur2 = [M c2 ; c2' z2];
    constraints = [constraints, Mschur1 >= 0, Mschur2 >= 0];
    objective = z1 + z2;

elseif (obj == "mrd+zep")

    c1 = [0; -1/(2*theta(3)); theta(2)/(2*theta(3)^2)];
    c2 = [0; -1/theta(3); theta(2)/theta(3)^2];
    z1 = sdpvar(1);
    z2 = sdpvar(1);
    Mschur1 = [M c1 ; c1' z1];
    Mschur2 = [M c2 ; c2' z2];
    constraints = [constraints, Mschur1 >= 0, Mschur2 >= 0];
    objective = z1 + z2;

elseif (obj == "rip")
    A = sqrt(theta(2)^2/(4*theta(3)^2) + log(1-p)/theta(3));
    % paper has mistake in gradient
    %g2 = theta(2)/(4*A) - 1/(2*theta(3));
    g2 = theta(2)/(2*theta(3)^2 * sqrt(theta(2)^2/theta(3)^2 + 4*log(1-p)/theta(3))) - 1/(2*theta(3));
    B = -theta(2)^2/(2*theta(3)^3) - log(1-p)/theta(3)^2;
    g3 = B/(2*A) + theta(2)/(2*theta(3)^2);
    c = [0; g2; g3];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;

elseif (obj == "rimp")
    g3 = theta(2)/(2*theta(3)^2) - log(1-p)/(2*theta(3)^2 * sqrt(log(1-p)/theta(3)));
    c = [0; -1/(2*theta(3)); g3];
    z = sdpvar(1);
    Mschur = [M c ; c' z];
    constraints = [constraints, Mschur >= 0];
    objective = z;

elseif (obj == "rip+rimp")
    % multiobjective
    % compute c-vector for rip
    A = sqrt(theta(2)^2/(4*theta(3)^2) + log(1-p)/theta(3));
    % paper has mistake in gradient
    %g2 = theta(2)/(4*A) - 1/(2*theta(3));
    g2 = theta(2)/(2*theta(3)^2 * sqrt(theta(2)^2/theta(3)^2 + 4*log(1-p)/theta(3))) - 1/(2*theta(3));
    B = -theta(2)^2/(2*theta(3)^3) - log(1-p)/theta(3)^2;
    g3 = B/(2*A) + theta(2)/(2*theta(3)^2);
    c1 = [0; g2; g3];

    % compute c-vector for rimp
    g3 = theta(2)/(2*theta(3)^2) - log(1-p)/(2*theta(3)^2 * sqrt(log(1-p)/theta(3)));
    c2 = [0; -1/(2*theta(3)); g3];

    % epigraph + Schur complement for both
    z1 = sdpvar(1);
    z2 = sdpvar(1);
    Mschur1 = [M c1 ; c1' z1];
    Mschur2 = [M c2 ; c2' z2];
    constraints = [constraints, Mschur1 >= 0, Mschur2 >= 0];

    % compound criterion
    objective = (1-alpha)*z1 + alpha*z2;

end


% call solver
optimize(constraints, objective, options);

% display results in a more readable format
if (type == "exact")
    design = [value(n) dose'];
    w = design(design(:, 1) > 0, 1);
    x = design(design(:, 1) > 0, 2);
elseif (type == "approx")
    design = [value(w) dose'];
    w = design(design(:, 1) >= 0.009, 1);
    x = design(design(:, 1) >= 0.009, 2);
end


end
