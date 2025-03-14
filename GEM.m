function [pert] = GEM(surf1, surf2, corr_true, options)

%% initialize 2-hop
n1 = surf1.n;
n2 = surf2.n;
EYE1 = sparse(1:n1, 1:n1, 1, n1, n1);
EYE2 = sparse(1:n2, 1:n2, 1, n2, n2);
W21 = sparse((double((surf1.adj * surf1.adj) > 0) - surf1.adj - EYE1) > 0);
W22 = sparse((double((surf2.adj * surf2.adj) > 0) - surf2.adj - EYE2) > 0);

%% initialize shots
pert = surf2.shots * surf1.shots';
pert = greedy_match(pert);

%% initialize filter
[~, U, lambda] = laplacian_from_TRIV_adj(surf2.adj, options.spec_dim);
[~, V, mu] = laplacian_from_TRIV_adj(surf1.adj, options.spec_dim);
Coeff = 1 ./ ((lambda - mu').^2 + 1);

%% local map improvement
[num, ~] = find(pert);
num = size(num, 1);
MA = full(diag(surf1.Phi.A));
e = zeros(num, 1);

%% Inside GEM function

% Check if 'surf2' has the 'distances' field
if ~isfield(surf2, 'distances')
    disp('Calculating geodesic distances for surf2...');
    
    % Define idx (starting point or set of points)
    idx = 1;  % or any other index based on your specific requirement
    
    % Ensure that idx is valid before passing to the fast marching function
    if any(idx > size(surf2.VERT, 1))
        error('idx is out of bounds for surf2.VERT');
    end

    % Compute distances using perform_fast_marching_mesh
    surf2.distances = perform_fast_marching_mesh(surf2.VERT, surf2.TRIV, idx, options.option2);
end

% Now you can safely access 'distances'
R_max = max(max(surf2.distances));

% Continue with the rest of the GEM function...


for kk = 1:options.maxIter
    pertF = greedy_match(W22 * pert * W21);
    [pertF, col] = find(pertF);
    DD1 = cell(surf1.n, 1);
    DD2 = cell(surf2.n, 1);
    idx = cell(num, 1);

    
    good = 1:min(size(surf1.distances, 1), size(surf2.distances, 1)); % Ensures both match


    corr_true = corr_true(corr_true <= size(surf1.distances, 1));

    % Correctly initialize 'ee' to be of length 'good'
    ee = zeros(length(good), 1);

    if any(good > size(surf2.distances, 1))
        error('ERROR: Index out of bounds for surf2.distances!');
    end





    % Safe indexing of surf1.distances and surf2.distances with corrected 'good' and 'pertF(good)'
    D1T = surf1.distances(corr_true(good), good);
    D2T = surf2.distances(good, good);  % Ensure it only gets needed rows/cols

  

    for i = 1:length(good)
        idx{i} = find(D1T(:, i) ~= 0);
        DD1{i} = D1T(idx{i}, i);
        valid_idx = idx{i}(idx{i} <= size(D2T, 1)); % Ensure indices are within bounds
        DD2{i} = D2T(valid_idx);
        r = max(DD1{i});
        ee(i) = sum(((abs(DD1{i} - DD2{i})) / r) .* MA(idx{i})) / sum(MA(idx{i}));
    end

    % Ensure 'ee' and 'good' have the same length before assignment
    if length(ee) ~= length(good)
        error('Mismatch in lengths of ee and good');
    end

    % Initialize 'e' if it has not been initialized already
    if kk == 1
        e = zeros(num, 1);  % Initialize 'e' to match the overall number of elements
    end

    e(good) = ee;  % Perform the assignment if lengths match

    if kk <= length(options.th)
        landmarks = find(e < options.th(kk));
    else
        landmarks = find(e < options.th(end));
    end
    sub_landmarks = landmarks;

    good = setdiff(1:num, sub_landmarks);
    goodF = setdiff(1:num, pertF(sub_landmarks));
    p = U(pertF(sub_landmarks), :)' * V(sub_landmarks, :);
    p = Coeff .* p;
    p = U(goodF, :) * p * V(good, :)';
    [p, ~] = greedy_match(p);
    [p, ~] = find(p);
    pertF(good) = goodF(p);
    pert = sparse(pertF, col, 1, surf2.n, surf1.n);
end
