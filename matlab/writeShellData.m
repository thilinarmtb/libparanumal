%% This code requires Chebfun.  To run it,
%%   1.  Download Chebfun:  git clone https://github.com/chebfun/chebfun
%%       or download and unpack the release *.zip/tarball.
%%   2.  In MATLAB, run addpath('/path/to/chebfun') to put Chebfun on the path.
%%   3.  Run this code by calling writeASBFDirichletData(N).
%%   4.  Put any asbfNXX.dat files generated in solvers/asbf/data/.

function writeShellData(N)
  fname = sprintf('shellN%02d.dat', N);
  fid = fopen(fname, 'w');

  x = chebfun(@(x) x);

  %%%%%%%%%% GLOBAL DISCRETE EIGENFUNCTION BASIS %%%%%%%%%%

  % Nodes and quadrature weights.
  [Rquad, Wquad] = legpts(N + 5);
  Wquad = Wquad.';
  Rgll = [-1 ; roots(diff(legpoly(N))) ; 1];
  Nplot = ceil(2.5*N);
  Rplot = linspace(-1, 1, Nplot)';

  writeFloatMatrix(fid, Rquad, 'SHELL GLOBAL DISCRETE QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad, 'SHELL GLOBAL DISCRETE QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rgll,  'SHELL GLOBAL DISCRETE GLL NODES');
  writeFloatMatrix(fid, Rplot, 'SHELL GLOBAL DISCRETE PLOT NODES');

  % BCs:  Dirichlet-Dirichlet
  T = (1 + x).*(1 - x).*chebpoly(0:(N - 1));
  DT = diff(T);

  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, VTquad,  'SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll,   'SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll,  'SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot,  'SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'SHELL GLOBAL DISCRETE INITIAL BASIS DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE');

  % BCs:  Neumann-Neumann
  T = [chebpoly(0) cumsum((1 + x).*(1 - x).*chebpoly(0:(N - 2)))];
  DT = diff(T);

  VTquad = feval(T, Rquad);
  VDTquad = feval(DT, Rquad);
  VTgll = feval(T, Rgll);
  VDTgll = feval(DT, Rgll);
  VTplot = feval(T, Rplot);
  VDTplot = feval(DT, Rplot);

  writeFloatMatrix(fid, VTquad,  'SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDTquad, 'SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTgll,   'SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN GLL VANDERMONDE');
  writeFloatMatrix(fid, VDTgll,  'SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, VTplot,  'SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDTplot, 'SHELL GLOBAL DISCRETE INITIAL BASIS NEUMANN-NEUMANN PLOT DERIVATIVE VANDERMONDE');

  %%%%%%%%%% PIECEWISE DISCRETE EIGENFUNCTION BASIS %%%%%%%%%%

  Rquad = [-1 ; roots(diff(legpoly(N + 4))) ; 1];
  T = chebpoly(0:(N + 4));
  V = T(Rquad);
  L = T*inv(V);
  Wquad = sum(L).';

  Nplot = ceil(2.5*N);
  Rplot = linspace(-1, 1, Nplot)';
  
  Rgll = [-1 ; roots(diff(legpoly(N))) ; 1];
  T = chebpoly(0:N);
  V = T(Rgll);
  L = T*inv(V);
  
  V = feval(L, Rquad);
  VD = feval(diff(L), Rquad);
  Vplot = feval(L, Rplot);
  VDplot = feval(diff(L), Rplot);

  writeFloatMatrix(fid, Rquad,  'SHELL PIECEWISE DISCRETE QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad,  'SHELL PIECEWISE DISCRETE QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rquad,  'SHELL PIECEWISE DISCRETE PLOT NODES');
  writeFloatMatrix(fid, Rgll,   'SHELL PIECEWISE DISCRETE GLL NODES');
  writeFloatMatrix(fid, V,      'SHELL PIECEWISE DISCRETE QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VD,     'SHELL PIECEWISE DISCRETE QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, Vplot,  'SHELL PIECEWISE DISCRETE PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDplot, 'SHELL PIECEWISE DISCRETE PLOT DERIVATIVE VANDERMONDE');
  
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

  writeFloatMatrix(fid, R,      'SHELL TRUE R');
  writeFloatMatrix(fid, lambda, 'SHELL TRUE LAMBDA');
  writeFloatMatrix(fid, Rquad,  'SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE NODES');
  writeFloatMatrix(fid, Wquad,  'SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE WEIGHTS');
  writeFloatMatrix(fid, Rgll,   'SHELL TRUE DIRICHLET-DIRICHLET GLL NODES');
  writeFloatMatrix(fid, Rplot,  'SHELL TRUE DIRICHLET-DIRICHLET PLOT NODES');
  writeFloatMatrix(fid, D,      'SHELL TRUE DIRICHLET-DIRICHLET EIGENVALUES');
  writeFloatMatrix(fid, Vquad,  'SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE VANDERMONDE');
  writeFloatMatrix(fid, VDquad, 'SHELL TRUE DIRICHLET-DIRICHLET QUADRATURE DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, Vgll,   'SHELL TRUE DIRICHLET-DIRICHLET GLL VANDERMONDE');
  writeFloatMatrix(fid, VDgll,  'SHELL TRUE DIRICHLET-DIRICHLET GLL DERIVATIVE VANDERMONDE');
  writeFloatMatrix(fid, Vplot,  'SHELL TRUE DIRICHLET-DIRICHLET PLOT VANDERMONDE');
  writeFloatMatrix(fid, VDplot, 'SHELL TRUE DIRICHLET-DIRICHLET PLOT DERIVATIVE VANDERMONDE');

  % TODO:  Neumann BCs.

  fclose(fid);
end
