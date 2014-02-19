function [x,relres,iter,resvec,mvecs] = mpgmres(A,b,P,type_in,tol,maxits,x0,varargin)
%MPGMRES An extension of GMRES which accepts multiple preconditioners.
%
%   X = MPGMRES(A,B,P) solves the linear system A*X = B for X.  The matrix
%   A should be of dimension n by n, and generally large and sparse.  The
%   vector B must have length n.  P should be a cell array of
%   preconditioners, i.e. n by n matrices (or function handles which act on
%   n-vectors) which approximate 'A' in some sense.  If the supplied P is a
%   single matrix, then the method is equivalent to right-preconditioned
%   GMRES.  
%
%   X = MPGMRES(A,B,P,TYPE) specifies the type of mpgmres algorithm used.
%   TYPE can be a string, in which case the accepted values are:
%    - 'full',which runs full MPGMRES, or 
%    - 'trunc', which runs a truncated version.
%   TYPE can also be a structure, which gives more control over the type of
%   truncation.  For details, see <a href="matlab: help mpgmres_120112>extended_help">here</a>.  
%
%   X = MPGMRES(A,B,P,TYPE,TOL) specifies the tolerance of the method.  The
%   default value is 1e-6.
%
%   X = MPGMRES(A,B,P,TYPE,TOL,MAXITS) specifies the maxiumum number of
%   iterations the methods allows.  The default value is n.
%
%   X = MPGMRES(A,B,P,TYPE,TOL,MAXITS,X0) specifies the initial guess.  X0
%   must be a vector of dimension n.  The default value is the zero vector.
%
%   X = MPGMRES(A,B,P,TYPE,TOL,MAXITS,X0,STOREZ) specifies whether to store
%   or calculate the matrix of search directions.  STOREZ = 1 stores, 
%   STOREZ = 0 calculates.  The default is 0. 
%
%   X = MPGMRES(A,B,P,TYPE,TOL,MAXITS,X0,STOREZ,TESTORTH) specifies whether
%   or not to check for loss of orthogonality (should only be used for
%   diagonostic purposes, as this is a relatively expensive process).
%   TESTORTH = 1 runs the check, TESTORTH = 0 does not.  The default is 0.
%
%   X = MPGMRES(A,B,P,TYPE,TOL,MAXITS,X0,STOREZ,TESTORTH,SAVEMATS)
%   specifies whether or not to save the matrices generated in the Arnoldi process at each iteration.  This
%   should only be used for diagonostic purposes, as it is a relatively 
%   expensive process.  SAVEMATS = 0 doesn't save matrices, SAVEMATS = 1
%   does.  The default is 0.
%
%   [X,RELRES] = MPGMRES(A,B,...) also returns the relative residual.
%
%   [X,RELRES,ITER] = MPGMRES(A,B,...) also returns the number of
%   iterations needed for convergence.
% 
%   [X,RELRES,ITER,RESVEC] = MPGMRES(A,B,...) also returns a vector of the
%   residual norms at each iteration.
%
%   [X,RELRES,ITER,RESVEC,MVECS] = MPGMRES(A,B,...) also returns a vector
%   of the cumulative number of inner products taken at the end of
%   each iteration 
%
%   For more information, see:
%   Greif, C., Rees, T., Szyld, D. B.,
%   <a href="http://www.cs.ubc.ca/~tyronere/TR-2011-12.pdf">Multi-preconditioned GMRES</a>
%   Technical report: UBC CS TR-2011-12, or Temple Math. report 11-12-23
%   
%   12 Jan 2012: First release version of MPGMRES
%
%   Tyrone Rees                        tyronere@cs.ubc.ca
%   Department of Computer Science 
%   University of British Columbia


[lb,wb] = size(b);
if (lb==1)||(wb==1) % test if  b is a vector
    if (wb>=1)&&(lb==1) % make sure it's a column vector
        b = b'; lb = wb;
    end
else
    error('right hand side b must be a vector')
end

if isa(A,'float')
    [lA,wA] = size(A);
    if lA ~= wA % test if A is a square matrix
        error('The input matrix A must be square')
    end    
    if lb ~= lA % test if the size of the matrix and the size of b
        error('The right hand side vector must be of the same length as A')
    end
    Atype = 'matrix';
elseif isa(A,'function_handle')
    Atype = 'func';
else
    error('unsupported type of A supplied')
end

if nargin < 3
    error ('Not enough input parameters: ensure A, b and P are specified')
end

if nargin < 4
    type_in = 'trunc';
end

if nargin < 5
    tol = 1e-6;
else
    if isa(tol,'numeric') ~= 1
        error('The tolerance must be a numeric value')
    end
end
if nargin < 6
    maxits = lb;
else
    if isa(maxits,'numeric') ~= 1
        error('The max number of iterations must be an integer')
    end
end
if nargin < 7
    x0 = zeros(lb,1);
else
    [lx0,wx0] = size(x0);
    if (lx0 == 1)&&(wx0 == lb)
        x0 = x0'; lx0 = wx0;
    end
    if lx0 ~= lb
        error('The starting vector must be of the same length as b')
    end
end
numvarargs = length(varargin);
if numvarargs > 4
    error('too many inputs')
end
optargs = {0 0 0}; % optional values of 'storeZ','testorth', 'savemats'
optargs(1:numvarargs) = varargin;
[storeZ, testorth, savemats] = optargs{:};
% testorth -- tests orthogonality, see Paige slide 14: 
%     http://www.stanford.edu/dept/ICME/docs/seminars/Paige-2009-11-04.pdf


%% Get the type and number of preconditoiners
if isa(P,'cell')
    k=length(P);  % number of preconditioners
else % otherwise, just one preconditioner, applies standard right pre-GMRES
    k = 1;  P = {P};
end

% check that it's a valid type and calculate the number of 'inner' interations for the maxit supplied
if isa(type_in,'char') % if just a string, then just a vanilla full or trunc mpgmres
    switch type_in
        case 'trunc' 
            nmaxits = maxits*k;
            type.type = 'trunc'; % set defaults
            type.col = 0;
            type.method = 'sum';
        case 'full' % full algorithm
            if k == 1
                nmaxits = maxits;
            else
                nmaxits = (maxits^(k+1) - maxits)/(maxits - 1);
            end
            type = type_in;
        otherwise
            error('Invalid type parameter')
    end
elseif isa(type_in,'struct') % if a structure, then allow more control
    type = type_in;
    switch type.type
        case 'trunc'
            nmaxits = maxits*k+1;
            switch type.col
                case 1
                    valid = {'inorder','reverse','alternate','random'};
                    if max(strcmp(type.method,valid)) == 0
                            error('invalid type.method parameter')
                    end
                case 0
                    valid = {'sum','random'};
                    if max(strcmp(type.method,valid)) == 0
                            error('invalid type.method parameter')
                    end
                otherwise
                    error('invalid type.col value: must be 0 or 1')
            end
    end
end

%% pre-allocate variables

indvar = min(lb,nmaxits);
if storeZ
    Z(lb,indvar) = 0.0
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %!! We don't need to store Z, as Z_k = [P_1^{-1} V_k ... P_t^{-1} V_k]%
    %!! and hence (Z^)y_k = P_1^{-1}(V^(:,Vindex_1))y_k(yindex_1) + ...   %
    %!!                     P_t^{-1}(V^(:,Vindex_t))y_k(yindex_t)         %
    %!! the values for the indices are stored in                          %
    %!! V_index and yk_index                                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize arrays
    end_V = cell(k,1); end_V(:) = {0};      % stores the end of Vk_index
    % initalize the structure that holds this information:
    Zinfo = struct('yk_index',zeros(nmaxits,1),'V_index',zeros(nmaxits,1),'end_V',end_V);
    lastV = 0;
end
V(lb,indvar+1) = 0.0;
H(indvar+1,indvar) = 0.0;
resvec(indvar,1) = 0.0;
mvecs(indvar,1) = 0.0;
c(indvar,1) = 0.0;
s(indvar,1) = 0.0;
rhs(indvar,1) = 0.0;

%% Set some flags
lindep_flag = 0;        % set the linear dependence flag - set to 1 if H_(p+1,p) = 0


%% Initialization
if strcmp(Atype,'matrix')
    r = b - A*x0;           % initial residual
else
    r = b - A(x0);
end
mvecs(1) = 1; mvecs(2) = mvecs(1);
nr = norm(r);           % norm of initial residual
rhs(1) = nr;            % initialize rhs
resvec(1) = 1;          % initialize the vector which stores the residuals
V(:,1) = r/nr;          % initialize the V matrix with normalized init. resid.
Z_temp = fullmultiprecondition(P,V(:,1));    % send to multiprecondition to get Z_temp
z_it = [1,0];           % [(outer) iteration, col of Z]
Zk = k;                 % a parameter which stores the size of the current Z_i
if ~storeZ
    VforZ = 1;
    ind_V = ones(k,1);
end
pp = 0;                 % initialize a concurrent index
nVk = 0;                % initialize the current size of Vk

%% Start the loop...
for p = 1:nmaxits                 % loop over the columns of Zk
    pp = pp+1;                    % update the 'real' index -- allows removal of rows
    if strcmp(Atype,'matrix')     % get vector which needs orthogonalizing
        w = A*Z_temp(:,z_it(2)+1);           % initial residual
    else
        w = A(Z_temp(:,z_it(2)+1));
    end      
    mvecs(z_it(1)+1) = mvecs(z_it(1)+1) + 1;
    for i = 1:pp                  % loop over columns of V
        H(i,pp) = V(:,i)'*w;
        w = w - H(i,pp)*V(:,i);   % remove component of V_i from w
    end
    nw = norm(w);
    H(pp+1,pp) = nw;          % calculate norm of vector w
    
    % TEST if w is the zero vector, hence Z_k is linearly dependent
    if H(pp+1,pp)<sqrt(eps); lindep_flag = 1;  end
    for l = 1:pp-1  % update previous rots
        h1 = H(l,pp);
        h2 = H(l+1,pp);
        H(l,pp) = c(l)*h1 + s(l)*h2;
        H(l+1,pp) = -s(l)*h1 + c(l)*h2;
    end
    h1 = H(pp,pp);
    h2 = H(pp+1,pp);
    gam = sqrt(h1^2 + h2^2);
    c(pp) = h1/gam;
    s(pp) = h2/gam;
    H(pp,pp) = c(pp)*h1+s(pp)*h2;
    H(pp+1,pp) = 0;
    
    % update rhs
    rhs_old = rhs;              % save old rhs
    rhs(pp+1) = -s(pp)*rhs(pp);
    rhs(pp) = c(pp)*rhs(pp);
    
    %    test_lindep_block;
    if (lindep_flag==1)&&((abs(rhs(pp+1)) >= tol*nr)||(isnan(abs(rhs(pp+1))))); % the last column of Z is dependent on the others
        fprintf('Column of Z linearly dependent...removing \n')
        pp = pp-1;          % reduce index by 1
        lindep_flag = 0;    % reset linear dependence flag
        rhs = rhs_old;      % replace rhs
    else
        nVk = nVk + 1;
        if storeZ
            Z(:,pp) = Z_temp(:,z_it(2)+1);
        else
            ik = ceil(nVk/VforZ); % find which preconditioner
            Zinfo(ik).end_V = Zinfo(ik).end_V + 1;
            Zinfo(ik).yk_index(Zinfo(ik).end_V) = pp;
            if isa(Zinfo(ik).V_index,'cell') % then type.col == 1, and cells are involved
                Zinfo(ik).V_index{Zinfo(ik).end_V} = ind_V;
                Zinfo(ik).wt{Zinfo(ik).end_V} = wt;
                lastV = max(lastV,max(Zinfo(ik).V_index{Zinfo(ik).end_V}));
            else
                Zinfo(ik).V_index(Zinfo(ik).end_V) = ind_V(nVk);
                lastV = max(lastV,Zinfo(ik).V_index(Zinfo(ik).end_V));
            end
            
        end
        if lindep_flag==1  % we have a lucky breakdown
            z_it(1) = z_it(1)+1;
            resvec(z_it(1)) = abs(rhs(pp+1))/nr;    % save residual
                                                    % to resvec
            fprintf('MPGMRES converged -- lucky breakdown! \n')
            break
        end
        V(:,pp+1) = w/nw;       % set new basis vector V_{p+1}
        switch testorth
            case 1             % test if V's are still orthogonal (see note above)
                S = triu(V'*V,1); fprintf('Measure of loss of orthog. \n')
                norm((eye(size(S)) + S)\S)
        end
    end
    z_it(2) = z_it(2)+1;         % update the index column of Z we're working on
    
    if z_it(2) == Zk
        resvec(z_it(1)+1) = abs(rhs(pp+1))/nr;    % save residual
                                                  % to resvec
        
        
        %% test convergence
        if resvec(z_it(1)+1)<tol
            fprintf('MPGMRES converged! \n')
            z_it(1) = z_it(1) + 1;
            break
        end
        %% multiprecondition
        if nVk == 0
            error('All current search directions are linearly dependent on the previous ones.\n Re-reun with a different trucation rule or starting vector.')
        end
        if isa(type,'char') % if just a string, then just a vanilla full or trunc mpgmres
            switch type
                case 'full' % full algorithm
                    ind_V = repmat(pp+1-nVk+1:pp+1,1,k);
                    Z_temp = fullmultiprecondition(P,V(:,pp+1-nVk+1:pp+1));
                    Zk = k*nVk;
                    if ~storeZ
                        VforZ = nVk;
                    end
                    nVk = 0; % reset nVk
            end
        elseif isa(type,'struct') % if a structure, then allow more control
            if type.col == 1 % act on individual columns of V
                % the supported methods in this case are currently:
                % in order, reverse order, random order, alternate order
                switch type.method
                    case 'inorder' % P_i^{-1}Ve_i
                        ind_V = pp+1 + (1-nVk:0);
                        Z_temp = multipreconditionvec(P,V(:,ind_V));
                        if ~storeZ
                            VforZ = 1;
                            ind_V = increase_ind_V(ind_V,k);
                        end
                    case 'reverse' % P_i^{-1}Ve_j, j = n-i+1
                        ind_V = fliplr(pp+1 + (1-nVk:0));
                        Z_temp = multipreconditionvec(P,V(:,ind_V));
                        if ~storeZ
                            VforZ = 1;
                            ind_V = increase_ind_V(ind_V,k);
                        end
                    case 'alternate' % swaps between previous two cases
                        if rem(pp,2) == 0
                            ind_V = fliplr(pp+1-nVk+1:pp+1);
                            Z_temp = multipreconditionvec(P,V(:,ind_V));
                            if ~storeZ
                                ind_V = increase_ind_V(ind_V,k);
                            end
                        else
                            ind_V = pp+1-nVk+1:pp+1;
                            Z_temp = multipreconditionvec(P,V(:,ind_V));
                            if ~storeZ
                                VforZ = 1;
                                ind_V = increase_ind_V(ind_V,k);
                            end
                        end
                    case 'random' % P_i^{-1}Ve_j, j random
                        order = randperm(nVk);
                        ind_V = zeros(1,nVk);
                        for iv = 1:nVk
                            ind_V(iv) = pp+1-nVk+order(iv);
                        end
                        Z_temp = multipreconditionvec(P,V(:,ind_V));
                        if ~storeZ
                            VforZ = 1;
                            ind_V = increase_ind_V(ind_V,k);
                        end
                    otherwise
                        error('Unsupported value of type.method')
                end
                Zk = k; nVk = 0;      % update size of Z_k and current size of V_k
            elseif type.col == 0 % act on all the columns of V
                if (z_it(1) == 1)&&(~storeZ) % upadate V_index so it's a cell array, not a vector                    
                    for ii = 1:k
                        Zinfo(ii).V_index = cell(nmaxits,1);
                        Zinfo(ii).V_index(1) = {1};
                        Zinfo(ii).wt = cell(nmaxits,1);
                        Zinfo(ii).wt(1) = {1};
                    end
                end 
                switch type.method     
                    case 'sum' % sum the columns
                        VforZ = 1;
                        ind_V = pp+1-nVk+1:pp+1;
                        wt = ones(nVk,1);
                        v_temp = V(:,ind_V)*wt;
                        Z_temp = fullmultiprecondition(P,v_temp);
                    case 'random' % sum the columns with random weights
                        ind_V = pp+1-nVk+1:pp+1;
                        wt = rand(nVk,1);
                        v_temp = V(:,ind_V)*wt;
                        Z_temp = fullmultiprecondition(P,v_temp);
                    otherwise
                        error('Unsupported value of type.method')
                end
                Zk = k; nVk = 0;      % update size of Z_k and current size of V_k
            else
                error('Please set type.col to be 0 or 1')
            end
        else
            error('Invalid ''type'' option. \n')
        end
        z_it(1) = z_it(1)+1;
        mvecs(z_it(1)+1)= mvecs(z_it(1));
        z_it(2) = 0;
        
    end
    switch savemats
        case 1
            save('mpgmres_const.mat','V','H','Z');
    end
end
switch savemats
    case 1
        save('mpgmres_b_const.mat','V','H','Z');
end

yk = H(1:pp,1:pp)\rhs(1:pp);
if storeZ
    x = x0 + Z(:,1:pp)*yk;
else
    x = x0;
    if isa(Zinfo(1).V_index,'cell') % then type.col = 0
        for ii = 1:k
            Vyk_i = zeros(lb,1);
            for iv = 1:Zinfo(ii).end_V % get the summed matrix
                Vyk_i = Vyk_i + (V(:,Zinfo(ii).V_index{iv})*Zinfo(ii).wt{iv})*yk(Zinfo(ii).yk_index(iv));
            end
            if isa(P{ii},'function_handle')
                x = x + P{ii}(Vyk_i);
            else
                x = x + P{ii}\Vyk_i;
            end
            
        end
    else
        for ii = 1:k
            Vyk_i = V(:,Zinfo(ii).V_index(1:Zinfo(ii).end_V))*yk(Zinfo(ii).yk_index(1:Zinfo(ii).end_V));
            if isa(P{ii},'function_handle')
                x = x + P{ii}(Vyk_i);
            else
                x = x + P{ii}\Vyk_i;
            end
        end
    end
end



resvec = resvec(1:z_it(1));
mvecs = mvecs(1:z_it(1));
iter = length(resvec) - 1;
relres = resvec(end)/norm(b);


%% apply multiple preconditioners to the residual r

function [z]=multipreconditionvec(pre,Q)

if isa(pre{1},'cell')
    [z]=multipreconditionvec(pre{1},Q);
else
    [n,m] = size(Q);
    k = length(pre);
    z=zeros(n,k);
    for i=1:m
        if isa(pre{i},'function_handle')
            z(:,i)=pre{i}(Q(:,i));
        else
            z(:,i)=pre{i}\Q(:,i);
        end
    end
    for i = m+1 : length(pre);
        if isa(pre{i},'function_handle')
            z(:,i)=pre{i}(Q(:,mod(i,m)+1));
        else
            z(:,i)=pre{i}\Q(:,mod(i,m)+1);
        end
    end
end

function [z]=fullmultiprecondition(pre,Q)

if isa(pre{1},'cell')
    [z]=fullmultiprecondition(pre{1},Q);
else
    lpre = length(pre);
    [lA,lQ] = size(Q);
    z=zeros(lA,lQ*lpre);
    ind = 1;
    ivec = 1:lQ*lpre;
    for i=1:lpre
        if isa(pre{i},'function_handle')
            for j = 1:lQ
                z(:,ivec(ind))=pre{i}(Q(:,j)); ind = ind+1;
            end
        else
            z(:,ivec(ind:ind+lQ-1))=pre{i}\Q; ind = ind + lQ;
        end
    end
end

function new_ind_V = increase_ind_V(ind_V,k)
sizeInd_V = length(ind_V);
if sizeInd_V < k
    new_ind_V = zeros(1,k);
    new_ind_V(1:sizeInd_V) = ind_V;
    for jv = sizeInd_V + 1:k
        new_ind_V(jv) = ind_V(mod(jv,sizeInd_V)+1);
    end
else
    new_ind_V = ind_V;
end


function extended_help
%EXTENDED_HELP some technical details about the 'type' option
%   As well as accepting the strings 'full' and 'trunc', MPGMRES also
%   accepts a structure input.  The structure should have the fields
%   TYPE.col and TYPE.method  Accepted values in this case are:
%    - TYPE.col = 1; these act on individual columns of the matrix of basis 
%      vectors, V_k.  Supported values of TYPE.method in this case are:
%       + TYPE.method = 'inorder'.  Applies P_i to the ith column of V_k +
%       TYPE.method = 'reverse'.  Applies P_i to the (t+1-i)th column of
%          V_k, where V_k has t columns
%       + TYPE.method = 'alternate'.  Alternates between the previous two
%          methods between iterations
%       + TYPE.method = 'random'.  Applies P_i to a random permutation of
%          V_k
%     - TYPE.col = 0; these act all of the columns of V_k.  Supported
%       values of TYPE.method in this case are:
%        + TYPE.method = 'sum'.  Applies P_i to the sum of the columns of
%           V_k
%        + TYPE.mehtod = 'random'. Applies P_i to a linear combination of
%           the columns of V_k, with weights uniformly distributed in the
%           range [0,1]
%   For more information, see the manuscript 
%   <a
%   href="http://www.cs.ubc.ca/~tyronere/TR-2011-12.pdf">Multi-preconditioned GMRES</a>.
error('This is a placeholder function just for helptext');

