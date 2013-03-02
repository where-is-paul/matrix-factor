%
%
% Main for mixed low-frequency Maxwell
% 
% we approximate solution (u,p) of
%
% curl curl u + nabla p = f, div u = 0     on the domain [-1,1]x[-1,1]  
% u x n = 0, p = 0  on boundary
%
% by FE-solution (uh, ph) in Ned(0,1) x P_1.
% Ned(0,1) denotes the FE-space of Nedelec edge-elements of first type and
% lowest degree.
% P_1 denotes the space of scalar nodal pw. linear elements.
%
% We create different types of plots for the computed FE-solution.
%
% Anna Schneebeli, June 02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% h = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.15];
ref=3;                      % number of global refinements
L2_error_u = zeros(1,ref);
Hcurl_error = zeros(1,ref);
L2_error_p = zeros(1,ref);
H1_error = zeros(1,ref);
l_num = zeros(1,ref);
DoFs = zeros(1,ref);
    
[p,e,t]=initmesh('squareg','hmax',1.5);
for i=1:ref
    
   clear te be AA BB LL A B L uh ph
    
    [p,e,t]=refinemesh('squareg',p,e,t,'regular');
    
    nod_num = length(p);          % number of nodes
  

    [te, be, e_num]=edges(t,e);
    [AA,BB,LL] = assemble_mixed(p,t,te,e_num);
    CC = zeros(nod_num, nod_num);
    [A,B,C,L] = boundary_cond(AA,BB,CC,LL,be,e);
    S=assemble_laplacian(p,t);
[S_bnd] = boundary_cond_laplacian(S,e);
    
  %h2=1/length(t);
  %h=sqrt(h2);
%       disp('Press any key for S_bnd') 
     %vec1(i)=condest(S_bnd);
       gamma1=normest(A)/normest(B)^2;
       gamma2=length(t);
       gamma3=1/length(t);
	gamma4=1 ;
    %vcond1(i)=condest(A+gamma1*B*B');
    vcond2(i)=cond(full(A+gamma2*B*B'))
    %vcond3(i)=condest(A+gamma3*B*B');
    %vcond4(i)=condest(A+gamma4*B*B');
    
    %vmin1(i)=eigs(A+gamma1*B*B',1,'SM');
    %vmin2(i)=eigs(A+gamma2*B*B',1,'SM');
    %vmin3(i)=eigs(A+gamma3*B*B',1,'SM');
    %vmin4(i)=eigs(A+gamma4*B*B',1,'SM');
    %vmax1(i)=eigs(A+gamma1*B*B',1,'LM');
    %vmax2(i)=eigs(A+gamma2*B*B',1,'LM');
    %vmax3(i)=eigs(A+gamma3*B*B',1,'LM');
    %vmax4(i)=eigs(A+gamma4*B*B',1,'LM');
    %if i>1 
	%i%,
	%vcond1(i)/vcond1(i-1),
	%%vcond2(i)/vcond2(i-1),
	%vcond3(i)/vcond3(i-1),
	%vcond4(i)/vcond4(i-1),
	%%vmin1(i)/vmin1(i-1),
	%vmin2(i)/vmin2(i-1),
	%vmin3(i)/vmin3(i-1),
	%vmin4(i)/vmin4(i-1),
	%vmax1(i)/vmax1(i-1),
	%vmax2(i)/vmax2(i-1),
	%vmax3(i)/vmax3(i-1),
	%vmax4(i)/vmin4(i-1),
    
    %end


% 
%  
%   [S_bnd] = boundary_cond_laplacian(S,e)
%    % solve linear system 
%    S = [A B;B' C];
%    F = [L;zeros(nod_num,1)];
%    sol = S\F;
%norm(sol),pause
%    uh = sol(1:e_num,1);
%    ph = sol(e_num+1:length(sol), 1);
%
%    [L2_error_u(i),Hcurl_error(i)] = uh_errors(p,t,te,uh);
%    [L2_error_p(i),H1_error(i)] = ph_errors(p,t,ph);
%
%    % compute number of elements
%    l_num(i) = length(t);
%    
%    % compute number of degrees of freedom in u
%    DoFs(i) = length(uh);
%
%    disp(['mesh ' num2str(i) ': done!'])
end


break;

figure(2)
    if i>1 i,vec1(i)/vec1(i-1),vec2(i)/vec2(i-1),end


% 
%  
%   [S_bnd] = boundary_cond_laplacian(S,e)
%    % solve linear system 
%    S = [A B;B' C];
%    F = [L;zeros(nod_num,1)];
%    sol = S\F;
%norm(sol),pause
%    uh = sol(1:e_num,1);
%    ph = sol(e_num+1:length(sol), 1);
%
%    [L2_error_u(i),Hcurl_error(i)] = uh_errors(p,t,te,uh);
%    [L2_error_p(i),H1_error(i)] = ph_errors(p,t,ph);
%
%    % compute number of elements
%    l_num(i) = length(t);
%    
%    % compute number of degrees of freedom in u
%    DoFs(i) = length(uh);
%
%    disp(['mesh ' num2str(i) ': done!'])
end


break;

figure(2)
    if i>1 i,vec1(i)/vec1(i-1),vec2(i)/vec2(i-1),end


% 
%  
%   [S_bnd] = boundary_cond_laplacian(S,e)
%    % solve linear system 
%    S = [A B;B' C];
%    F = [L;zeros(nod_num,1)];
%    sol = S\F;
%norm(sol),pause
%    uh = sol(1:e_num,1);
%    ph = sol(e_num+1:length(sol), 1);
%
%    [L2_error_u(i),Hcurl_error(i)] = uh_errors(p,t,te,uh);
%    [L2_error_p(i),H1_error(i)] = ph_errors(p,t,ph);
%
%    % compute number of elements
%    l_num(i) = length(t);
%    
%    % compute number of degrees of freedom in u
%    DoFs(i) = length(uh);
%
%    disp(['mesh ' num2str(i) ': done!'])
end


break;

figure(2)
subplot(2,1,1)
loglog(l_num,L2_error_u,'o-',l_num,Hcurl_error,'*-',l_num,l_num.^(-1/2),'r')
% loglog(l_num,Hcurl_error,'*-',l_num,l_num.^(-1/2),'r')
title('L2- and H(curl)-error of uh vs. number of elements')
legend('L2-error','H(curl)-error','Order 1')
%title('H(curl)-error vs. number of elements')
%legend('H(curl)-error','Order 1')
xlabel('number of elements in mesh')
ylabel('norm-errors')
grid on

subplot(2,1,2)
loglog(l_num,L2_error_p,'o-',l_num,H1_error,'*-',l_num,l_num.^(-1/2),l_num,l_num.^(-1),'-.' )
% loglog(l_num,Hcurl_error,'*-',l_num,l_num.^(-1/2),'r')
title('L2- and H1-error of ph vs. number of elements')
legend('L2-error','H1-error','Order 1', 'Order 2')
%title('H(curl)-error vs. number of elements')
%legend('H(curl)-error','Order 1')
xlabel('number of elements in mesh')
ylabel('norm-errors')
grid on



% plots of the solution (uh, ph)
%
x=linspace(-1,1,15);
y=x;

% surface plot of uh
u_fem = uh_on_grid(p,t,te,uh,x,y);
disp('here we go! 2-D plots...')
z1=zeros(length(y),length(x));
z2=z1;
z1(:,:)=u_fem(1,:,:);
z2(:,:)=u_fem(2,:,:);

% plot of x-component of uh
figure(3)
% pdemesh(p,e,t,zeros(1,length(p)))
% hold on
mesh(x,y,z1)
title('x-component of FE-solution')
xlabel('x')
ylabel('y')
hidden on

% plot of y-component of uh
figure(4)
surf(x,y,z2)
title('y-component of FE-solution')
xlabel('x')
ylabel('y')


% vector-field plot of uh
disp('here we go! vector plots...')
figure(5)
quiver(x,y,z1,z2)
title('Vector-field plot of FE-solution uh')
xlabel('x')
ylabel('y')

% plot of the solution ph on finest grid with pde-toolbox commands
figure(6)
pdeplot(p,e,t, 'xydata', ph, 'zdata', ph, 'mesh', 'on')
title('FE-pressure ph')


% Here we automatically save the data in a file called [filename].mat
% The data can be extracted by typing the command >> load [filename]  in the MATLAB command window
%save example2

%return
