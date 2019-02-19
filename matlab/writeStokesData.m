% TODO:  Remove Chebfun dependency.

function writeStokesData(N)
	fname = sprintf('stokes%02d.dat', N);
	fid = fopen(fname, 'w');

	Q = legpoly(0:N, 'norm');

	xk = [-1 ; roots(diff(Q(:, end))) ; 1];
	yk = [-1 ; roots(diff(Q(:, end - 1))) ; 1];

	V = Q(xk);
	I = eye(N + 1, N + 1);
	I(N + 1, N + 1) = 0.0;

	W = Q(yk);

	P1 = V*I*inv(V);  
	P2 = feval(Q(:, 1:end-1), xk)*inv(feval(Q(:, 1:end-1), yk))*feval(Q, yk)*inv(V);
	P3 = feval(Q(:, 1:end-1), xk)*pinv(feval(Q(:, 1:end-1), xk));

	writeFloatMatrix(fid, P1, 'Pressure projection matrix - L2-orthogonal');
	writeFloatMatrix(fid, P2, 'Pressure projection matrix - Interpolatory');
	writeFloatMatrix(fid, P3, 'Pressure projection matrix - Pseudoinverse');

	fclose(fid);
end
