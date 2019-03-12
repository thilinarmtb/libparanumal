% This code requires Chebfun.
%
% TODO:  Remove this dependency.
function writeStokesData(N)
	fname = sprintf('stokesN%02d.dat', N);
	fid = fopen(fname, 'w');

	cubNq = ceil(3*N/2) + 1;
	[cubxk, cubwk] = legpts(cubNq);
	cubwk = cubwk.';
	cubQ = legpoly(0:(cubNq - 1));
	cubV = feval(cubQ, cubxk);
	cubL = cubQ*inv(cubV);
	cubD = feval(diff(cubL), cubxk);

	writeFloatMatrix(fid, cubD, 'Cubature 1D differentiation matrix');

	NP = N - 1;
	xkP = [-1 ; roots(diff(legpoly(NP))) ; 1];
	QP = legpoly(0:NP);
	VP = QP(xkP);
	LP = QP*inv(VP);
	cubInterpP = LP(cubxk);

	writeFloatMatrix(fid, cubInterpP, 'Pressure 1D cubature interpolation matrix');

	fclose(fid);
end
