function c = c_rauc(theta)
% compute gradient of ratio of areas under curve function

    % compute c-vector for RAUC objective
    % evaluate RAUC for grid of theta values
    theta1_grid = [theta(1), linspace(theta(1)-0.1,theta(1)+0.1, 20)];
    theta2_grid = [theta(2), linspace(theta(2)-0.1,theta(2)+0.1, 20)];
    theta3_grid = [theta(3), linspace(theta(3)-0.1,theta(3)+0.1, 20)];
    [theta1, theta2, theta3] = ndgrid(theta1_grid, theta2_grid, theta3_grid);
    points = [theta1(:), theta2(:), theta3(:)];
    values = arrayfun(@(i) rauc(points(i,:)), 1:size(points,1));
    F_grid = reshape(values, size(theta1));

    % compute gradient
    [g1, g2, g3] = gradient(F_grid);

    % extract value of gradient at point
    t = (theta1 == theta(1)) & (theta2 == theta(2)) & (theta3 == theta(3));
    indt = find(t);
    c = [g1(indt); g2(indt); g3(indt)];
    
end