function density = ssdMatrixGenerator(dim, maxAbsCoeff, definite, sign, filename)
% SSDMATRIXGENERATOR  Create a symmetric (semi-) definite matrix
%
% This function creates and writes positve or negative
% (semi-)definite matrices and is used to set up the QP data for
% some of the numerical experiments in the paper "Computing
% Feasible Points of Bilevel Problems with a Penalty Alternating
% Direction Method" by Thomas Kleinert and Martin Schmidt.
%
% Syntax: density = ssdMatrixGenerator(dim, maxAbsCoeff, definite, sign, filename)
%
% Inputs:
%    dim         - dimension of the output matrix
%    maxAbsCoeff - maximum absolute value of the matrix entries
%    definite    - 0 = semidefinite, 1 = strictly definite
%    sign        - -1 = negative definite, 1 = positive definite
%    filename    - filename of the output matrix file
%
% Outputs:
%    density - density of the ouput matrix
%
%
% Author(s): Thomas Kleinert, Martin Schmidt
% Email: thomas.kleinert@fau.de, martin.schmidt@uni-trier.de
% Website: https://github.com/m-schmidt-math-opt/qp-bilevel-matrix-generator


%------------- BEGIN CODE --------------------------------------------------------
tic
rng('shuffle');

%---------------------------------------------------------------------------------
% handle inputs
assert(dim >= 1);
assert(definite == 0 || definite == 1);
assert(sign == 1 || sign == -1);
if maxAbsCoeff == 0
  maxAbsCoeff = 1;
end

%---------------------------------------------------------------------------------
% 1. generate random matrix S with size dim x dim and density baseDensity
%    i) compute baseDensity
dims = [0, 25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, ...
        15000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, ...
        125000, 150000, 175000, 200000, 250000, 500000, 1000000, 1500000, ...
        2000000];
densities = [1.0, 1.0, 0.5, 0.2, 0.15, 0.06, 0.036, 0.027, 0.022, 0.0135, ...
             0.009, 0.0072, 0.006, 0.0045, 0.0035, 0.0025, 0.0015, ...
             0.0011, 0.0008, 0.0006, 0.0005, 0.0004, 0.00025, 0.00022, ...
             0.00018, 0.00014, 0.0001, 0.00005, 0.00001, 0.000005, ...
             0.000001, 0.000001];
i = 0;
for d = dims
    if d >= dim
        break;
    else
        i = i + 1;
    end
end
slope = (densities(i) - densities(i+1))/(dims(i) - dims(i+1));
baseDensity = densities(i) + slope*(dim - dims(i));
%    ii) generate vectors rows, cols, and values
nonzeros = int32(baseDensity*dim*dim);
rows = randi(dim, nonzeros, 1);
cols = randi(dim, nonzeros, 1);
values = randi([-maxAbsCoeff,maxAbsCoeff], nonzeros, 1);
%    iii) replace 0-entries in values by random integer
values(values==0) = randi(maxAbsCoeff);
%    iv) make sure that there is an entry for the last row and column
posr = randi(nonzeros);
rows(posr,1) = dim;
posc = randi(nonzeros);
cols(posc,1) = dim;
%    v) build sparse matrix S
S = sparse(rows, cols, values);

%---------------------------------------------------------------------------------
% 2. generate a positive (semi-)definite matrix M from S
%    i) compute M = S'*S
M = S'*S;
%    ii) make strictly definite if parameter is set accordingly
if definite == 1
    d = randi([1, maxAbsCoeff], dim, 1);
    M = M + spdiags(d(:),0,dim,dim);
end
%    iii) Set correct sign
M=sign*M;

%---------------------------------------------------------------------------------
% 3. write M to file
%    i) get the triplets
[row, col, v] = find(M);
%    ii) write out the matrix
data_dump = [row, col, v];
fid = fopen(filename, 'w');
fprintf( fid, '%d %d %d\n', transpose(data_dump));
fclose(fid);

%---------------------------------------------------------------------------------
% compute actual density to return
density = nnz(M)/(dim*dim);
fprintf(" > Wrote matrix with density %f to file.\n", density);

toc
end
%------------- END OF CODE -------------------------------------------------------
