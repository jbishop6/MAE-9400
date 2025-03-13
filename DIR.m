function pertF = DIR(name1, name2, options, corr_true)

addpath(genpath('.'));
spec_dim = options.spec_dim; % max spectrum dimension
spec_dim_cut = options.spec_dim_cut; 
option1.nb_iter_max = 30; % iteration number needed for fast marching
option2.nb_iter_max = 120;
th = options.th;
iter_number = options.maxIter;

%%
% [surf2.pt, surf2.trg] = ReadOFF(name2);
data2 = load(strcat(name2, '.mat'));  
if isfield(data2, 'fv') && isfield(data2.fv, 'vertices') && isfield(data2.fv, 'faces')
    surf2.pt = data2.fv.vertices;
    surf2.trg = data2.fv.faces;
else
    error('.mat file for surf2 does not contain expected structure.');
end

S2 = MESH.MESH_IO.read_shape(name2);
surf2.Phi = MESH.compute_LaplacianBasis(S2, spec_dim);
surf2.n = length(surf2.pt);
num2 = surf2.n;
opts.shot_num_bins = 10; % number of bins for shot
opts.shot_radius = 5; % percentage of the diameter used for shot
Xdesc = calc_shot(surf2.pt', surf2.trg', 1:num2, opts.shot_num_bins, opts.shot_radius * surf2.Phi.sqrt_area / 100, 3)';
vertex2 = surf2.pt';
faces2 = surf2.trg';

% [surf1.pt, surf1.trg] = ReadOFF(name1);
data1 = load(strcat(name1, '.mat'));  
if isfield(data1, 'fv') && isfield(data1.fv, 'vertices') && isfield(data1.fv, 'faces')
    surf1.pt = data1.fv.vertices;
    surf1.trg = data1.fv.faces;
else
    error('.mat file for surf1 does not contain expected structure.');
end

S1 = MESH.MESH_IO.read_shape(name1);

surf1.Phi = MESH.compute_LaplacianBasis(S1, spec_dim);  
surf1.n = length(surf1.pt);
MA = full(diag(surf1.Phi.A));
xdesc = calc_shot(surf1.pt', surf1.trg', 1:surf1.n, opts.shot_num_bins, opts.shot_radius * surf1.Phi.sqrt_area / 100, 3)';
pertF = knnsearch(Xdesc, xdesc, 'NSMethod', 'kdtree');

if surf1.n < surf2.n
    num = surf1.n;
else
    num = surf2.n;
end
vertex1 = surf1.pt';
faces1 = surf1.trg';
cnt = 0;
e = zeros(num, 1);
D1 = perform_fast_marching_mesh(vertex1, faces1, 1:surf1.n, option1);
D2 = perform_fast_marching_mesh(vertex2, faces2, 1:surf2.n, option2);
R_max = max(max(D2));
landmarks = [];

for kk = 1:iter_number
    good = 1:num;  % Initialize 'good' inside the loop
    valid_good = good(good <= size(D1, 2));  % Ensure 'good' is within bounds of D1

    % Debugging print statements
    disp('Iter:');
    disp(kk);
    disp('Size of D1:');
    disp(size(D1));
    disp('Size of D2:');
    disp(size(D2));
    disp('Values in corr_true(good):');
    disp(corr_true(valid_good));
    disp('Values in good:');
    disp(valid_good);
    disp('Length of valid_good:');
    disp(length(valid_good));
    disp('Length of valid_corr_true:');
    disp(length(corr_true(valid_good)));

    % Correctly initialize 'ee' based on 'valid_good'
    ee = zeros(length(valid_good), 1);

    % Safe indexing of D1 and D2 with corrected 'good' and 'pertF(good)'
    D1T = D1(corr_true(valid_good), valid_good);
    D2T = D2(valid_good);  % Corrected indexing here

    disp('D1T size:');
    disp(size(D1T));
    disp('D2T size:');
    disp(size(D2T));

    DD1 = cell(length(valid_good), 1);
    DD2 = cell(length(valid_good), 1);
    idx = cell(length(valid_good), 1);

    for i = 1:length(valid_good)
        idx{i} = find(D1T(:, i) ~= 0);
        DD1{i} = D1T(idx{i}, i);
        DD2{i} = D2T(idx{i});  % Correctly indexing 'D2T' as a column vector
        
        % Handle empty arrays
        if isempty(DD1{i}) || isempty(DD2{i})
            disp('Empty array detected: skipping this iteration.');
            ee(i) = NaN;  % Assign NaN to skip this iteration
            continue;
        end
        
        DD2{i}(DD2{i} == 0) = R_max;

        % Print the sizes and values of DD1{i} and DD2{i} to debug dimension mismatch
        disp('Size of DD1{i}:');
        disp(size(DD1{i}));
        disp('Size of DD2{i}:');
        disp(size(DD2{i}));
        disp('Values in DD1{i}:');
        disp(DD1{i});
        disp('Values in DD2{i}:');
        disp(DD2{i});

        % Check if sizes match before performing the calculation
        if length(DD1{i}) ~= length(DD2{i})
            error('Mismatch in dimensions of DD1{i} and DD2{i}');
        end

        r = max(DD1{i});
        if r == 0
            r = 1; % Avoid division by zero
        end

        % Calculate ee(i) safely with matching dimensions
        ee(i) = sum(((abs(DD1{i} - DD2{i})) / r) .* MA(idx{i})) / sum(MA(idx{i}));
    end

    % Remove NaN values from ee before assignment
    ee = ee(~isnan(ee));
    valid_good = valid_good(~isnan(ee));

    % Debugging lengths before assignment
    disp('Length of ee:');
    disp(length(ee));
    disp('Length of valid_good:');
    disp(length(valid_good));

    % Ensure 'ee' and 'valid_good' have the same length before assignment
    if length(ee) ~= length(valid_good)
        error('Mismatch in lengths of ee and valid_good');
    end

    e(valid_good) = ee;  % Perform the assignment if lengths match

    if kk <= length(th)
        landmarks = find(e < th(kk));
    else
        landmarks = find(e < th(end));
    end
    sub_landmarks = landmarks;
    good = setdiff(1:num, sub_landmarks);
    goodF = setdiff(1:num, pertF(sub_landmarks));
    H = surf1.Phi.evecs(sub_landmarks, :)' * surf2.Phi.evecs(pertF(sub_landmarks), :);
    [U, D, V] = svd(H);
    dim = findK(diag(D));
    if dim > spec_dim_cut
        cnt = cnt + 1;
    end
    if cnt > 3
        break;
    end
    specB = surf2.Phi.evecs(:, 1:dim);
    specA = surf1.Phi.evecs(:, 1:dim);
    H = specA(sub_landmarks, :)' * specB(pertF(sub_landmarks), :);
    [U, D, V] = svd(H);
    C = U * V';
    specBB = specB(goodF, :);
    specAA = specA(good, :);
    p = knnsearch(specBB * C', specAA, 'NSMethod', 'kdtree');
    pertF(good) = goodF(p);
end
