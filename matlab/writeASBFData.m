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

  %%%%%%%%%% GLOBAL DISCRETE EIGENFUNCTION BASIS %%%%%%%%%%

  % Nodes and quadrature weights.
  [Rquad, Wquad] = legpts(N + 5);
  Wquad = Wquad.';
  Rgll = [-1 ; roots(diff(legpoly(N))) ; 1];
  Nplot = ceil(2.5*N);
  Rplot = linspace(-1, 1, Nplot)';

  writeFloatMatrix(fid, Rquad, 'ASBF GLOBAL DISCRETE QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad, 'ASBF GLOBAL DISCRETE QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rgll,  'ASBF GLOBAL DISCRETE GLL NODES');
  writeFloatMatrix(fid, Rplot, 'ASBF GLOBAL DISCRETE PLOT NODES');

  % BCs:  Dirichlet-Dirichlet
  T = (1 + x).*(1 - x).*chebpoly(0:(N - 1));
  DT = diff(T);

  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, VTquad,  'ASBF GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'ASBF GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll,   'ASBF GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll,  'ASBF GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot,  'ASBF GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'ASBF GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE');

  % BCs:  Neumann-Neumann
  T = [chebpoly(0) cumsum((1 + x).*(1 - x).*chebpoly(0:(N - 2)))];
  DT = diff(T);

  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, VTquad,  'ASBF GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'ASBF GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll,   'ASBF GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll,  'ASBF GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot,  'ASBF GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'ASBF GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN PLOT DERIVATIVE VANDERMONDE');

  %%%%%%%%%% TRUE EIGENFUNCTION BASIS %%%%%%%%%%

  % TODO:  Is it possible to remove the hard-coded dependence on R and lambda?
  % (Need to understand how R and lambda affect the true eigenvalues.)
  R      = 1.5;
  lambda = 1.0;

  % BCs:  Dirichlet-Dirichlet
  L = chebop(1, R);
  L.op = @(r, u) r.^2.*diff(u, 2) + 2*r.*diff(u) - lambda*r.^2.*u;
  L.lbc = 0;
  L.rbc = 0;

  [Q, D] = eigs(L, N);
  Q = quasimatrix(fliplr(Q));
  D = -flipud(diag(D));

  for (i = 1:1:N)
    Q(:, i) = Q(:, i)/norm(Q(:, i));
  end

  Nquad = length(Q) + 2;
  [Rquad, Wquad] = legpts(Nquad, [1 R]);
  Wquad = Wquad.';
  Rgll = (1 + [-1 ; roots(diff(legpoly(N))) ; 1])*((R - 1)/2) + 1;
  Nplot = ceil(2.5*length(Q));
  Rplot = linspace(1, R, Nplot)';

  DQ = diff(Q);

  Vquad = feval(Q, Rquad);
  VDquad = feval(DQ, Rquad);
  Vgll = feval(Q, Rgll);
  VDgll = feval(DQ, Rgll);
  Vplot = feval(Q, Rplot);
  VDplot = feval(DQ, Rplot);

  writeFloatMatrix(fid, R,      'ASBF TRUE R');
  writeFloatMatrix(fid, lambda, 'ASBF TRUE LAMBDA');
  writeFloatMatrix(fid, Rquad,  'ASBF TRUE DIRICHLET-DIRICHLET QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad,  'ASBF TRUE DIRICHLET-DIRICHLET QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rgll,   'ASBF TRUE DIRICHLET-DIRICHLET GLL NODES');
  writeFloatMatrix(fid, Rplot,  'ASBF TRUE DIRICHLET-DIRICHLET PLOT NODES');
  writeFloatMatrix(fid, D,      'ASBF TRUE DIRICHLET-DIRICHLET EIGENVALUES');
  writeFloatMatrix(fid, Vquad,  'ASBF TRUE DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDquad, 'ASBF TRUE DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, Vgll,   'ASBF TRUE DIRICHLET-DIRICHLET GLL VANDERMONDE');
  writeFloatMatrix(fid, VDgll,  'ASBF TRUE DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, Vplot,  'ASBF TRUE DIRICHLET-DIRICHLET PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDplot, 'ASBF TRUE DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE');

  % TODO:  Neumann BCs.

  fclose(fid);
end
