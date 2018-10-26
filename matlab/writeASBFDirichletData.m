%% This code requires Chebfun.  To run it,
%%   1.  Download Chebfun:  git clone https://github.com/chebfun/chebfun
%%       or download and unpack the release *.zip/tarball.
%%   2.  In MATLAB, run addpath('/path/to/chebfun') to put Chebfun on the path.
%%   3.  Run this code by calling writeASBFDirichletData(N).
%%   4.  Put any asbfNXX.dat files generated in solvers/asbf/data/.

function writeASBFDirichletData(N)
  fname = sprintf('asbfN%02d.dat', N);
  fid = fopen(fname, 'w');

  % Generate a basis that satisfies the Dirichlet conditions.
  %
  % TODO:  Do the same for Neumann and write out the associated matrices.
  x = chebfun(@(x) x);
  T = (1 + x).*(1 - x).*chebpoly(0:(N - 1));
  DT = diff(T);

  [Rquad, Wquad] = legpts(N + 5);
  Wquad = Wquad.';
  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);

  Rgll = [-1 ; roots(diff(legpoly(N))) ; 1];
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);

  Nplot = ceil(2.5*N);
  Rplot = linspace(-1, 1, Nplot)';
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, Rquad, 'ASBF QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad, 'ASBF QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rgll, 'ASBF GLL NODES');
  writeFloatMatrix(fid, Rplot, 'ASBF PLOT NODES');

  writeFloatMatrix(fid, VTquad, 'ASBF INITIAL BASIS QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'ASBF INITIAL BASIS QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll, 'ASBF INITIAL BASIS GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll, 'ASBF INITIAL BASIS GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot, 'ASBF INITIAL BASIS PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'ASBF INITIAL BASIS PLOT DERIVATIVE VANDERMONDE');

  fclose(fid);
end
