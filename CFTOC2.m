function [xopt,uopt,flag] = CFTOC2(dyn,E,uf,vf,Q,R,x0,xf,xlb,xub,ulb,uub,dT,N)
disp('initialize...')

M=10^9;
%states:[x ;y ;v ;m ;theta]
%input:[T ;yaw]
x = sdpvar(5,N+1);
u = sdpvar(2,N);

%cost
cost=[];
%initial states constriant:
constraint = [x(:,1)==x0,x(1:3,end)==xf(1:3)];  %initial/final state

for k=1:N
    %states dynamics
    %LTV+wind:
    %constraint = constraint + [x(:,k+1)==Ak(x(:,k),u(:,k))*x(:,k)+Bk(x(:,k))*u(:,k)+E*[uf(x(1,k),x(2,k),k);vf(x(1,k),x(2,k),k)]];
    
    %LTV
    %constraint = constraint + [x(:,k+1)==Ak(x(:,k),u(:,k))*x(:,k)+Bk(x(:,k))*u(:,k)];
    constraint = constraint + [x(:,k+1)==dyn(x(:,k),u(:,k))];
    
    %LTI+wind
    %constraint = constraint + [x(:,k+1)==Ak*x(:,k)+Bk*u(:,k)+E*[uf(x(1,k),x(2,k),k);vf(x(1,k),x(2,k),k)]];
    
    %LTI
    %constraint = constraint + [x(:,k+1)==Ak*x(:,k)+Bk*u(:,k)];
    
    %states constraints
    constraint = constraint + [xlb <= x(:,k)<= xub];
    constraint = constraint + [-0.6 <= (x(3,k+1)-x(3,k))/dT<= 0.6];
    
    %input constraints
    constraint = constraint + [ulb <= u(:,k)<= uub];
    %cost
    cost=cost+(x(:,k)-xf)'*Q*(x(:,k)-xf)+u(:,k)'*R*u(:,k);
end

%cost = cost+ M*x(3,end);

%slove optimization problem
disp('optimizing...')
%choose solver
options=sdpsettings('solver','ipopt','verbose',0,'showprogress',1);
%options.ipopt.check_derivatives_for_naninf='yes';
result=optimize(constraint,cost,options);

%check feasibility
flag=result.problem;
disp(yalmiperror(flag));
if flag==0 
    %solution found
    xopt=double(x);    uopt=double(u);
else
    %no feasible solution
    xopt=[];    uopt=[];
end
end





