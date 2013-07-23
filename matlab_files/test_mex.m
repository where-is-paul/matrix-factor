%   Tests the factorization generated by matrix-factor mex file.
% 
%   No parameters are needed. Currently this script only tests
%   the following matrices:
%       aug3dcqp.mtx
%       bloweya.mtx
%       bratu3d.mtx
%       tuma1.mtx
%       tuma2.mtx
%       1138_bus.mtx
%
%   Outputs for each test case:
%       The fill factor (nnz(L+D+L')/nnz(A))
%       The number of iterations for GMRES to converge.

display('==============starting tests================');
warning off;
% testmats = [];
% for i = 1:8
%     testmats = [testmats; strcat('testmat', num2str(i))];
% end
% testmats = cellstr(testmats);
% 
% Lshape_mats = [];
% for i = 1:4
%     Lshape_mats = [Lshape_mats; strcat('Lshape_matrices', num2str(i))];
% end
% Lshape_mats = cellstr(Lshape_mats);
% 
% homogenous_mats = [];
% for i = 1:4
%     homogenous_mats = [homogenous_mats; strcat('matrices_homogeneous', num2str(i))];
% end
% homogenous_mats = cellstr(homogenous_mats);

other_mats = { 'aug3dcqp'; 'bloweya'; 'bratu3d'; ...
            'tuma1'; 'tuma2'; '1138_bus'; 'Lshape_matrices4'; 'lund_b'};

%other_mats = {'m3d10-001'};

%aug3dcqp has a terrible condition number
all_mats = other_mats;%[testmats; Lshape_mats; homogenous_mats; other_mats];
        
lfil = 2.2;
tol = 0.001;
for i = 1:length(all_mats)
    mat_name = all_mats{i};
    fprintf('Now testing %s:\n', mat_name);
    base = '';
    file = strcat(base, mat_name, '.mtx');
    A = mmread(file);
    
    [l d p S B] = ildl(A, lfil, tol);
   
    fprintf('The relative residual is %d.\n', norm(B - l*d*l', 1)/norm(B, 1));
    fprintf('The fill factor is %.3f.\n', nnz(l+d+l')/nnz(B));
    fprintf('The largest elem. of L is %.3f.\n', full(max(max(abs(l)))));
	fprintf('A has %i nnz.\n', nnz(A));
    %fprintf('The condition number is %d.\n', condest(B));

    %spy(A(p,p)); figure; spy(abs(l*d*l') > 1e-8);
    %spy(B); figure; spy(A(p,p));
    
    e = ones(size(B,1),1);
    %[y, flag, relres, iter, resvec] = ...
        gmres(B,S^(-1)*e,min(60,size(B,1)),1e-8,3,l*d, l');
    
    %semilogy(1:length(resvec), resvec, 'r-');

    %xlabel('iteration');
    %ylabel('relative residual');
    fprintf('\n');
end

warning on;
display('================ending tests================');
