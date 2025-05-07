function R = rauc(theta)
% calculate the ratio of areas under the curve
% used in computing gradient for c-optimality criterion
% numerically estimate AUC
m = @(x, theta1, theta2, theta3) exp(theta1 + theta2*x + theta3*x.^2);
I = integral(@(x) m(x, theta(1), theta(2), theta(3)), 0, -theta(2)/theta(3));

R = 1 + (exp(theta(1))*theta(2)/theta(3))/I;

end