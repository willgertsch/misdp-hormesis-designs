function M = compute_M(x,w,theta)
%compute_M Compute information matrix for design

    grad = @(x) [
        ones(1, length(x)); 
        x'; 
        x'.^2
        ]';
    eta = theta(1) + theta(2) .* x + theta(3) .* x.^2;
    v = exp(eta);
    F = grad(x);
    M = zeros(3,3);
    for i=1:length(x)
        M = M + w(i) * F(i, :)' * F(i, :) * v(i);
    end
end