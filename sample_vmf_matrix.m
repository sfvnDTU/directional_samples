function samples = sample_vmf_matrix(C, varargin)
%% Gibbs Sampler for the von Mises-Fischer Matrix distribution
% --- Description...??
%   Input: 
%
%   Output: 
%
% Based on following articles:
%       [1] Hoff, P. D. (n.d.). 
%           "Simulation of the matrix Bingham-von Mises-Fisher distribution, 
%           with applications to multivariate and relational data."
%           Retrieved from https://www.stat.washington.edu/research/reports/2007/tr526.pdf
%       [2] Wood, Andrew TA. "Simulation of the von Mises Fisher distribution." 
%           Communications in statistics-simulation and computation 23.1 (1994): 157-164.
% 
% Written by. Søren Føns Vind Nielsen, Nov. 2016, sfvn@dtu.dk

% Get inputs
opts = mgetopt(varargin);
n_samples = mgetopt(opts,'samples',100);
thinning = mgetopt(opts,'thinning',1); % save all samples


[p,N] = size(C);
assert(p>N) % sampler as of now only works in case of more observations than dimensions

[U,S,V] = svd(C,'econ');
X = U*V';

maxiter = thinning*n_samples;
samples = nan(p,N, n_samples);
ns = 1;

for i = 1:maxiter
    % Sample each column individually - random order
    columns = randperm(N);
    for c = columns
        cidx = true(1,N);
        cidx(c) = false;
        Nspace =  null(X(:,cidx)); % basis for nullspace of X withouth c'th column
        z = sample_vmf_vector(Nspace*C(:,c));
        X(:,c) = Nspace*z;
    end
    
    % Save
    if mod(i,thinning)==0
       samples(:,:,ns) = X; 
    end
    
    % Print statistics
        % Eval vMF matrix distribution and display correlation between samples?
end


%eof
end

% -------------------------------------------------------------------------
function out = mgetopt(varargin)
% MGETOPT Parser for optional arguments
%
% Usage
%   Get a parameter structure from 'varargin'
%     opts = mgetopt(varargin);
%
%   Get and parse a parameter:
%     var = mgetopt(opts, varname, default);
%        opts:    parameter structure
%        varname: name of variable
%        default: default value if variable is not set
%
%     var = mgetopt(opts, varname, default, command, argument);
%        command, argument:
%          String in set:
%          'instrset', {'str1', 'str2', ... }
%
% Example
%    function y = myfun(x, varargin)
%    ...
%    opts = mgetopt(varargin);
%    parm1 = mgetopt(opts, 'parm1', 0)
%    ...

% Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

if nargin==1
    if isempty(varargin{1})
        out = struct;
    elseif isstruct(varargin{1})
        out = varargin{1}{:};
    elseif isstruct(varargin{1}{1})
        out = varargin{1}{1};
    else
        out = cell2struct(varargin{1}(2:2:end),varargin{1}(1:2:end),2);
    end
elseif nargin>=3
    opts = varargin{1};
    varname = varargin{2};
    default = varargin{3};
    validation = varargin(4:end);
    if isfield(opts, varname)
        out = opts.(varname);
    else
        out = default;
    end
    
    for narg = 1:2:length(validation)
        cmd = validation{narg};
        arg = validation{narg+1};
        switch cmd
            case 'instrset',
                if ~any(strcmp(arg, out))
                    fprintf(['Wrong argument %s = ''%s'' - ', ...
                        'Using default : %s = ''%s''\n'], ...
                        varname, out, varname, default);
                    out = default;
                end
            case 'dim'
                if ~all(size(out)==arg)
                    fprintf(['Wrong argument dimension: %s - ', ...
                        'Using default.\n'], ...
                        varname);
                    out = default;
                end
            otherwise,
                error('Wrong option: %s.', cmd);
        end
    end
end
%eof
end