function S = compute_adjacencyBasis(S,numEigs)
%Compute the eigen vectors of the adjacency matrix
if nargin < 2, numEigs = 200; end
fprintf('Computing %d adjacency_Eigen-functions...',numEigs); tic;
A = S.adj; %+ zeros(S.n);
A =sparse(double(A));
disp(size(A))
% try
[evecs, evals] = eigs(A, numEigs);
% catch
%     % In case of trouble make the adjacency definite but would need to
%     add self loop corresponding to smallest eigval
%     [evecs, evals] = eigs(A + 1e-8*speye(S.n), numEigs, 'sm');
% end
evals = diag(evals);
% if ~isreal(evecs) %% will be real since A is real and symmetric
%     evecs(1:2:end) = real(evecs(1:2:end));
%     evecs(2:2:end) = imag(evecs(1:2:end));
% end
% % [evals, order] = sort(abs(evals),'ascend');
S.evals = evals;
% % evecs = evecs(:,order);
S.evecs = evecs;
t = toc; fprintf('done:%.4fs\n',t);
end