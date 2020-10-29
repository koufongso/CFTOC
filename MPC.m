clear;
clc;
%% Load Data
load uf;
load vf;
tracks=csvread('VRD211.csv',1,1);

%%
ORD=[272.093, 41.9742]; %start coord
SFO=[237.617, 37.6211]; %final coord   

%% Simulation Parameters
dT=60;          %[s] time step
N=4*3600/dT;    %Total time step

%% Drag Force Parameters:
Cd= 0.04;               %Drag force coefficient[-]
rho= 0.38;              %0.379637582    %[kg/m3] air density
eta= 0.01667*10^-3;     %0.596+0.581   0.5885[kg/(kN*s)]...SFC -> [kg/(N*s)]
Area= 122.6;            %[m2]...Aircraft wing area

%% Wind Data
%actual wind
%ua=@(x,y) interp1(Lon,Lat,ua,x,y,'linear');
%va=@(x,y) interp1(Lon,Lat,va,x,y,'linear');

Lon=235:1:275;
Lat=45:-1:35;

%forcast wind
ufun=@(x,y,t) interp2(Lon,Lat,squeeze(uf(:,:,t)),x,y,'linear');
vfun=@(x,y,t) interp2(Lon,Lat,squeeze(vf(:,:,t)),x,y,'linear');

%% Dynamic model/Linearized model
%Aircraft dynamics
%x=[x;y;v;m;theta]      u=[T;yaw];
%{
dyn=@(x,u) [x(1)+(dT*x(3)*cos(x(5)))/(111189.3*cos(x(2)));
            x(2)+(dT*x(3)*sin(x(5)))/111189.3;
            x(3)+dT*(2*u(1)-Cd*rho*Area*(x(3)^2))/(2*x(4));
            x(4)-dT*eta*u(1);
            x(5)+dT*u(2)];
%}

dyn=@(x,u) [x(1)+(dT*x(3)*cos(x(5)))/(111189.3*cos(x(2)));
            x(2)+(dT*x(3)*sin(x(5)))/111189.3;
            x(3)+dT*(2*u(1)-Cd*rho*Area*(x(3)^2))/(2*x(4));
            x(4)-dT*eta*u(1);
            x(5)+dT*9.81*u(2)/(x(3)+0.1)];
        
        
%Linearized model
%x(k+1)=Ak*x(k)+Bk(xk)+E*w(k)

Ak=@(x,u) eye(5)+diag([dT/(111493*cos(x(2))),dT/111111,dT,dT,dT])*[0,0,cos(x(5)),0,-232.5*sin(x(5));
                                                                 0,0,sin(x(5)),0,232.5*cos(x(5));
                                                                 0,0,-(Cd*rho*Area*232.5)/x(4),-(2*u(1)-Cd*rho*Area*(232.5^2))/(2*(x(4)^2)),0;
                                                                 0,0,0,0,0;
                                                                 0,0,-9.81*u(2)/(x(3)^2),0,0];
         
Bk=@(x) diag([dT/(111493*cos(x(2))),dT/111000,dT,dT,dT])*[0,0;0,0;1/x(4),0;-eta,0;0,9.81/x(3)];

E=[1,0; 0,1; 0,0; 0,0; 0,0];

%LTI
xeq=[(ORD+SFO)'/2; 232.5; 57615; deg2rad(190)];
ueq=[109000*2;0];
A=Ak(xeq,ueq);
B=Bk(xeq);

%A=@(x,u) eye(5)+diag([dT/(111493*cos(x(2))),dT/111000,dT,dT,dT])*[0,0,deg2rad(190),0,-232.5*sin(deg2rad(190));
                                                                   %0,0,sin(deg2rad(190)),0,232.5*cos(deg2rad(190));
                                                                   %0,0,-(Cd*rho*Area*232.5)/57615,-(2*100000*2-Cd*rho*Area*(232.5^2))/(2*(57615^2)),0;
                                                                   %0,0,0,0,0;
                                                                   %0,0,0,0,0];

%B=@(x) diag([dT/(111493*cos(x(2))),dT/111000,dT,dT,dT])*[0,0;0,0;1/57615,0;-eta,0;0,1];

%% Cost Matrics  
Q=diag([1000,1000,10000,0,0]);    %state cost
R=diag([0,0]);              %input cost

%% Initial/Final States
%x=[x; y; v; m; theta]      u=[T; yaw];
x0=[ORD'; 0; 78000; deg2rad(180)];      %initial states
xf=[SFO'; 0; 0; 0];   %final states

%% State/Input Constraints
%xlb=[235; 35; 555/3.6; 37230; -inf];
%xub=[275; 45; 871/3.6; 78000; inf];
xlb=[235; 35; 0; 37200; -inf];
xub=[275; 45; 250; 78000; inf];
ulb=[0; -0.6];
uub=[120000*2; 0.6];

%% States/Input Record
%xopt=nan(5,N+1);
%uopt=nan(2,N);
%xsim=nan(2,Tf+1);

%xopt(:,1)=x0;
%xsim(:,1)=x0;

%% Path Planning
%CFTOC(Ak,Bk,E,uf,vf,Q,R,x0,xf,xlb,xub,ulb,uub,N)
%{
tic
[xopt,uopt,flag]=CFTOC2(dyn,E,ufun,vfun,Q,R,x0,xf,xlb,xub,ulb,uub,dT,N);
toc
%}

%% LTI Path Planning
%CFTOC(Ak,Bk,E,uf,vf,Q,R,x0,xf,xlb,xub,ulb,uub,N)

tic
[xopt,uopt,flag]=CFTOC(A,B,E,ufun,vfun,Q,R,x0,xf,xlb,xub,ulb,uub,N);
toc


%% MPC Controller & Path Following
%{
MPCt=120; %[s] when to run MPC
for k=1:Tf
    %call MPC 
    if ~mod(k,MPCt)
        [xtemp,utemp]=CFTOC(Ak(x(:,k),u(:,k)),Bk(x(:,k)),wuf,wvf,Q,R,x(:,k));
        uopt(:,k:k+MPCt-1)=utemp(:,1:MPCt);
        xopt(:,k:k+MPCt-1)=xtemp(:,1:MPCt);
        counter=0;
    end
    
    %simulation
    xsim(:,k)=dyn(x(:,k),uopt(:,k));
end
%}

%% Visualization
%trajectory
figure;
title('flight trajectory')
plot(xopt(1,:),xopt(2,:),'-*');
hold on
plot(tracks(:,2)+360,tracks(:,1),'-o');
xlabel('lon');  ylabel('lat')
legend('optimized trajectory','actual trajectory')
grid on

%states-longtitude
figure
subplot(5,1,1)
plot(0:N,xopt(1,:));
hold on
plot(0:size(tracks,1)-1,tracks(:,2)+360)
legend('optimized','actual')
xlabel('k');    ylabel('lon');
grid on

%states-latitude
subplot(5,1,2)
plot(0:N,xopt(2,:));
hold on
plot(0:size(tracks,1)-1,tracks(:,1))
legend('optimized','actual')
xlabel('k');    ylabel('lat');
grid on   

%states-speed
subplot(5,1,3)
plot(0:N,xopt(3,:));
hold on
plot(0:size(tracks,1)-1,tracks(:,3)*0.44704)
legend('optimized','actual')
xlabel('k');    ylabel('v[m/s]');
grid on    

%states-mass
subplot(5,1,4)
plot(0:N,xopt(4,:)/1000)
xlabel('k');    ylabel('m[ton]');
grid on    

%states-"heading angle"
subplot(5,1,5)
plot(0:N,rad2deg(xopt(5,:)))
xlabel('k');    ylabel('theta[deg]');
grid on    

%input
figure
subplot(2,1,1)
plot(0:N-1,uopt(1,:)/1000)
xlabel('k');    ylabel('T [kN]');
grid on  
 
subplot(2,1,2)
plot(0:N-1,uopt(2,:))
xlabel('k');    ylabel('yaw');
grid on  

%%





