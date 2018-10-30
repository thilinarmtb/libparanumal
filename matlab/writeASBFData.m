%% This code requires Chebfun.  To run it,
%%   1.  Download Chebfun:  git clone https://github.com/chebfun/chebfun
%%       or download and unpack the release *.zip/tarball.
%%   2.  In MATLAB, run addpath('/path/to/chebfun') to put Chebfun on the path.
%%   3.  Run this code by calling writeASBFDirichletData(N).
%%   4.  Put any asbfNXX.dat files generated in solvers/asbf/data/.

function writeASBFData(N)
  fname = sprintf('asbfN%02d.dat', N);
  fid = fopen(fname, 'w');

  x = chebfun(@(x) x);

  % Nodes and quadrature weights.
  [Rquad, Wquad] = legpts(N + 5);
  Wquad = Wquad.';
  Rgll = [-1 ; roots(diff(legpoly(N))) ; 1];
  Nplot = ceil(2.5*N);
  Rplot = linspace(-1, 1, Nplot)';

  writeFloatMatrix(fid, Rquad, 'ASBF QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad, 'ASBF QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rgll,  'ASBF GLL NODES');
  writeFloatMatrix(fid, Rplot, 'ASBF PLOT NODES');

  % BCs:  Dirichlet-Dirichlet
  T = (1 + x).*(1 - x).*chebpoly(0:(N - 1));
  DT = diff(T);

  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, VTquad,  'ASBF INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'ASBF INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll,   'ASBF INITIAL BASIS DIRICHLET-DIRICHLET GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll,  'ASBF INITIAL BASIS DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot,  'ASBF INITIAL BASIS DIRICHLET-DIRICHLET PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'ASBF INITIAL BASIS DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE');

  % BCs:  Neumann-Neumann
  T = [chebpoly(0) cumsum((1 + x).*(1 - x).*chebpoly(0:(N - 2)))];
  DT = diff(T);

  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, VTquad,  'ASBF INITIAL BASIS NEUMANN-NEUMANN QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'ASBF INITIAL BASIS NEUMANN-NEUMANN QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll,   'ASBF INITIAL BASIS NEUMANN-NEUMANN GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll,  'ASBF INITIAL BASIS NEUMANN-NEUMANN GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot,  'ASBF INITIAL BASIS NEUMANN-NEUMANN PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'ASBF INITIAL BASIS NEUMANN-NEUMANN PLOT DERIVATIVE VANDERMONDE');

  fclose(fid);
end
