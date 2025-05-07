% alpha_grid = 0.1:0.1:0.9
% p = 0.50
function eff_tab = RI_dual(alpha_grid, p)
% Find dual-objective designs for estimating RI_p and RIM_p
% leaving a few things hard coded for now just to produce results
    theta = [3.086, 4.844,-13.113];

    % compute designs for alpha=0,1
    [x0, w0] = find_design_hormesis_count(100, 1, theta, 23, "exact", ...
        "rip", p, 0, 360);
    [x1, w1] = find_design_hormesis_count(100, 1, theta, 23, "exact", ...
        "rimp", p, 1.0, 360);
    
    % compute information matrices for each
    M0 = compute_M(x0, w0, theta);
    M1 = compute_M(x1, w1, theta);

    % compute c-vectors
    A = sqrt(theta(2)^2/(4*theta(3)^2) + log(1-p)/theta(3));
    g2 = theta(2)/(2*theta(3)^2 * sqrt(theta(2)^2/theta(3)^2 + 4*log(1-p)/theta(3))) - 1/(2*theta(3));
    B = -theta(2)^2/(2*theta(3)^3) - log(1-p)/theta(3)^2;
    g3 = B/(2*A) + theta(2)/(2*theta(3)^2);
    c0 = [0; g2; g3];

    g3 = theta(2)/(2*theta(3)^2) - log(1-p)/(2*theta(3)^2 * sqrt(log(1-p)/theta(3)));
    c1 = [0; -1/(2*theta(3)); g3];

    % compute objective values
    obj0 = c0' * (M0 \ c0);
    obj1 = c1' * (M1 \ c1);

    % loop through grid of alpha values
    eff_tab = zeros(length(alpha_grid), 4);
    for i = 1:length(alpha_grid)

        % find design with given alpha
        [xi,wi] = find_design_hormesis_count(100, 1, theta, 23, "exact", ...
            "rip+rimp", p, alpha_grid(i), 360);
        
        % compute information matrix
        Mi = compute_M(xi,wi,theta);
        
        % compute efficiency
        rip_eff = obj0 / (c0'*(Mi\c0));
        rimp_eff = obj1 / (c1'*(Mi\c1));

        % store
        eff_tab(i, :) = [alpha_grid(i), rip_eff, rimp_eff, length(xi)];
    end

    % add rows for alpha=0,1
    eff_tab = [
        0, obj0/(c0'*(M0\c0)), obj1/(c1'*(M0\c1)), length(x0);
        eff_tab;
        1, obj0/(c0'*(M1\c0)), obj1/(c1'*(M1\c1)), length(x1)
        ];

end