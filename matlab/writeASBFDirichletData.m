% This code requires Chebfun.  To run it,
%   1.  Download Chebfun:  git clone https://github.com/chebfun/chebfun
%       or download and unpack the release *.zip/tarball.
%   2.  In MATLAB, run addpath('/path/to/chebfun') to put Chebfun on the path.
%   3.  Run this code by calling writeASBFDirichletData(N).
%   4.  Put any asbfNXX.dat files generated in solvers/asbf/data/.

function writeASBFDirichletData(N)
	R = 1.5;  % TODO:  This should probably be an input parameter.

	fname = sprintf('asbfN%02d.dat', N);
	fid = fopen(fname, 'w');

	% Discrete Dirichlet eigenfunctions.
	%
	% TODO:  Remove Chebfun dependency.
	r = chebfun(@(r) r, [1 R]);
	U = (1 - r).*(R - r).*legpoly(0:N, [1 R], 'norm');
	A = diff(U)'*((r.^2).*diff(U));
	B = U'*(r.^2.*U);
	[V, D] = eig(A, B, 'vector');
	[D, p] = sort(D);
	V = V(:, p);
	Q_DDir = U*V;
	D_DDir = D;

	% Normalize in the r^2-weighted inner product.
	for (i = 1:1:size(Q_DDir, 2))
		qi = Q_DDir(:, i);
		normqi = sqrt(qi.'*(r.^2.*qi));
		Q_DDir(:, i) = qi/normqi;
	end

	[Rquad, Wquad] = legpts(N + 5, [1 R]);
	Wquad = (Rquad.^2).*(Wquad.');
	Bquad = feval(Q_DDir, Rquad);

	lambda = D_DDir;
	writeFloatMatrix(fid, lambda, 'ASBF EIGENVALUES');

	writeFloatMatrix(fid, Bquad, 'ASBF QUADRATURE VANDERMONDE');
	writeFloatMatrix(fid, Rquad, 'ASBF QUADRATURE NODES');
	writeFloatMatrix(fid, Wquad, 'ASBF QUADRATURE WEIGHTS');

	Rgll = [1 ; roots(diff(legpoly(N, [1 R]))) ; R];
	Bgll = feval(Q_DDir, Rgll);

	writeFloatMatrix(fid, Bgll, 'ASBF GLL VANDERMONDE');
	writeFloatMatrix(fid, Rgll, 'ASBF GLL NODES');

	fclose(fid);
end
