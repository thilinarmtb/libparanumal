% This code requires Chebfun.
%
% TODO:  Remove this dependency.
function writeStokesData(N)
	fname = sprintf('stokesN%02d.dat', N);
	fid = fopen(fname, 'w');

	Q = legpoly(0:N, 'norm');

	xk = [-1 ; roots(diff(Q(:, end))) ; 1];
	yk = [-1 ; roots(diff(Q(:, end - 1))) ; 1];

	W = Q(xk);
	I = eye(N + 1, N + 1);
	I(N + 1, N + 1) = 0.0;

	P1 = W*I*inv(W);  
	[U, S, V] = svd(eye(N + 1) - P1);
	u1 = S(1, 1)*U(:, 1);
	v1 = V(:, 1);

	P2 = feval(Q(:, 1:end-1), xk)*inv(feval(Q(:, 1:end-1), yk))*feval(Q, yk)*inv(W);
	[U, S, V] = svd(eye(N + 1) - P2);
	u2 = S(1, 1)*U(:, 1);
	v2 = V(:, 1);

	P3 = feval(Q(:, 1:end-1), xk)*pinv(feval(Q(:, 1:end-1), xk));
	[U, S, V] = svd(eye(N + 1) - P3);
	u3 = S(1, 1)*U(:, 1);
	v3 = V(:, 1);

	writeFloatMatrix(fid, u1, 'Pressure projection rank-1 update - L2 orthogonal (u)');
	writeFloatMatrix(fid, v1, 'Pressure projection rank-1 update - L2 orthogonal (v)');
	writeFloatMatrix(fid, u2, 'Pressure projection rank-1 update - Interpolatory (u)');
	writeFloatMatrix(fid, v2, 'Pressure projection rank-1 update - Interpolatory (v)');
	writeFloatMatrix(fid, u3, 'Pressure projection rank-1 update - Pseudoinverse (u)');
	writeFloatMatrix(fid, v3, 'Pressure projection rank-1 update - Pseudoinverse (v)');

	cubNq = ceil(3*N/2) + 1;
	[cubxk, cubwk] = legpts(cubNq);
	cubwk = cubwk.';
	cubQ = legpoly(0:(cubNq - 1));
	cubV = feval(cubQ, cubxk);
	cubL = cubQ*inv(cubV);
	cubD = feval(diff(cubL), cubxk);

	writeFloatMatrix(fid, cubD, 'Cubature 1D differentiation matrix');

	fclose(fid);
end
