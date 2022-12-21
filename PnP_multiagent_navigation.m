close all
clear all
clc

%% Desired Parameters of Recorded Output
global WantMovie WantSplitFrames frameNum frameRate Nframes timeBtwFrames

WantMovie=1;
WantSplitFrames=0;      % Only set this to 1 if you want to save a sequence of separate frames to the disk
frameRate = 30;         % Frame rate for playing the video
Nframes = 900;          % Total number of frames in the video


%% Select simulation parameters

global independentNav softSpring gamma adjMat varrho mu omega alpha

independentNav=0; % Does everyone navigate independently? (1 only needed for a single demo illustrating independent motion to the target. Otherwise make sure this is set to 0)
softSpring=0;   % Set to 1 for lazy-PnP and to 0 for contractive-PnP 

alpha = 1;        % PnP Controller order
mu = 0.015;        % PnP Control gain (in the paper, this corresponds to $\mu$ for the contractive case, and to $\omega$ in the lazy case).
gamma = 1;       % Leader Control gain

flowTime = 300; % Define simulation flow time
jumpTime = 1;   % Define the total number of jumps
maxStep = .1;    % Maximum step size for numerical solver
timeBtwFrames = flowTime/Nframes;


%% Construct the Environment and Target

%global variables for workspace polytope description
global workspaceIsSphereWorld workspaceMatrix workspaceCoefficients interiorWorkspacePoint

%select type of workspace
workspaceIsSphereWorld=true;

% global variables for spherical obstacle description
global stateDim obstacleNum obstacleCenters obstacleRadii obstacleClearance obstacleNumPoints

% other global variables in environment
global agentNum leader xStar R
agentNum=15;
leader=1;
xStar=[14,14]'; % initialize target
R=2; % communication radius

stateDim=2; % Dimension of state-space
obstacleNum=3; % Number of obstacles
obstacleNumPoints=360; % Number of vertices to represent the obstacle in graphics display
obstacleCenters=zeros(obstacleNum,stateDim); % Centers of reference balls, one for each obstacle, TBD
obstacleRadii=zeros(obstacleNum,1); % Radii of reference balls, one for each obstacle, TBD
obstacleClearance=zeros(obstacleNum,1);
for i=1:obstacleNum
    % set [circular] obstacle centers, interior radii, and clearances 
    theta=2.*pi*(i-1)/obstacleNum;
    obstacleCenters(i,:)=8*[cos(theta) sin(theta)];
    obstacleRadii(i)=i; % obstacle internal radii for the star-convex case
end
if workspaceIsSphereWorld
    obstacleRadii=[5;4;3]; % only good for the sim with 3 disk-shaped obstacles!
end

% Set clearances for buffer regions around the star-shaped obstacles
% This is best done manually: adjust the clearances to make sure the buffers of different obstacles DO NOT intersect. 
obstacleClearance=9*[0.65; 0.65; 0.5]; % only good for 3 obstacles!

% Construct the workspace:
% The workspace is the rectangle [wXmin,wXmax]\times[wYmin,wYmax]
global wXmin wXmax wYmin wYmax
wXmax=15; wXmin=-15; wYmax=15; wYmin=-15;
workspaceMatrix=[1 0; -1 0; 0 1; 0 -1]; 
workspaceCoefficients=[wXmax; -wXmin; wYmax; -wYmin];
interiorWorkspacePoint=[0; 0];

% Convert workspace into polygon form (list of vertices)
[workspaceVertices,nonRedundantFaces,nonRedundantEqualities]=qlcon2vert(interiorWorkspacePoint,workspaceMatrix,workspaceCoefficients);
workspacePolygon=convhull(polyshape(workspaceVertices,'SolidBoundaryOrientation','ccw','Simplify',false,'KeepCollinearPoints',true));
obsDisks=polyshape();
for j=1:obstacleNum
    obstacle=polybuffer([obstacleCenters(j,:)],'points',obstacleRadii(j));
    obsDisks=union(obsDisks,obstacle);
end
if workspaceIsSphereWorld
    workspacePolygon=subtract(workspacePolygon,obsDisks);
else
    angles=[linspace(0,double(2*pi),obstacleNumPoints),0];
    obstacleBuffer=polyshape();
    for j=1:obstacleNum
        %rad=barrierCurve(j);
        oc=obstacleCenters(j,:);
        radii=barrierCurve(j,angles);
        %compute obstacle as polygon
        [x,y]=pol2cart(angles,radii);
        x=x+oc(1);
        y=y+oc(2);
        obstacle=polyshape(x,y,'SolidBoundaryOrientation','ccw','Simplify',false,'KeepCollinearPoints',true);
        workspacePolygon=subtract(workspacePolygon,obstacle);
        %compute obstacle buffer as polygon
        radii=sqrt(barrierCurve(j,angles).^2+obstacleClearance(j)^2);
        [x,y]=pol2cart(angles,radii);
        x=x+oc(1);
        y=y+oc(2);
        obstacle=intersect(workspacePolygon,polyshape(x,y,'SolidBoundaryOrientation','ccw','Simplify',false,'KeepCollinearPoints',true));
        %obstacle=scale(obstacle,obstacleClearance(j)/obstacleRadii(j)+1,obstacleCenters(j,:));
        obstacleBuffer=union(obstacleBuffer,obstacle);
    end
end

%% Initialize agent positions
x0 = zeros(agentNum,stateDim);
%divide the agents into three "spans" of roughly equal size
dtheta=pi/(agentNum);
rems=mod(agentNum,3);
quot=(agentNum-rems)/3;
theta=-5*pi/6+quot*dtheta;
for k = 0:quot-1
    theta=theta-dtheta;
    x0(k+1,:)=obstacleCenters(1,:)+(obstacleRadii(1)+3)*[cos(theta),sin(theta)];
end
theta=pi+theta;
for k = quot:2*quot+rems-1
    theta=theta+dtheta;
    x0(k+1,:)=obstacleCenters(3,:)+(obstacleRadii(3)+2)*[cos(theta),sin(theta)];
end
theta=2*pi-theta;
for k=2*quot+rems:agentNum-1
    theta=theta-dtheta;
    x0(k+1,:)=obstacleCenters(2,:)+(obstacleRadii(2)+5)*[cos(theta),sin(theta)];
end

% Reshape the initial condition for the hybrid systems solver
x0Alt = reshape(x0',[1 agentNum*stateDim]);

%% Define the Adjacency matrix

% Generate an undirected adjacency matrix for the communication graph and compute the maximum distance along an edge
distMat = zeros(agentNum,agentNum);
R_max = 0;
for i = 1:agentNum
    for j = 1:agentNum
        if (norm(x0(i,:)-x0(j,:)) <= R) && (i ~= j)
            distMat(i,j) = 1;
        end
        
        if (norm(x0(i,:)-x0(j,:)) > R_max) && (distMat(i,j)==1)
           R_max =  norm(x0(i,:)-x0(j,:));
        end
    end
end

circMat = zeros(agentNum,agentNum);
for i = 1:agentNum
    for j = 1:agentNum
        if abs(i-j) == 1 || abs(i-j) == agentNum-1
            circMat(i,j) = 1;
        else
            circMat(i,j) = 0;
        end
    end
end

intMat = circMat;
broken = 0;      % Select which edge to remove from the circle
if broken <= agentNum-1 && broken >= 1
   intMat(broken,broken+1) = 0;
   intMat(broken+1,broken) = 0;
else
    intMat(1,agentNum) = 0;
    intMat(agentNum,1) = 0;
end

adjMat = intMat; % <-------------------------------------- Pick your adjacency matrix

EE = sum(adjMat*ones(agentNum,1))/2; % Number of edges

%Compute a value for the buffer radius
InitMaxDistBelowR=0;
for i = 1:agentNum
    for j = 1:agentNum
        if distMat(i,j)==1 && adjMat(i,j)==1
            InitMaxDistBelowR=max(InitMaxDistBelowR,norm(x0(i,:)-x0(j,:)));
        end
    end
end
varrho = 0.5*(InitMaxDistBelowR+R); % computed buffer radius                             

%% set control gains
m = R / varrho;
if softSpring==0
    varrhostar=varrho;
else
    varrhostar=varrho+(R-varrho)/(EE^(1/(2.0+alpha)));
end
M = r(R) / r(varrhostar);

% Define omega (edge tension function parameter)
if softSpring==0
    % This is the contractive case. The parameter omega is in charge of
    % guaranteeing the condition for graph maintenance, as opposed to the
    % parameter mu, which merely determines the gain associated with the
    % PnP field
    omega=mu*((2+alpha)/(2*R^(1+alpha)))*(m^(1+alpha)*(EE-m^2))/((m-1)^(2+alpha));
else
    % This is the Lazy PnP case. The parameter omega no longer influences
    % the graph maintenance guarantee and serves merely as a gain on the
    % PnP controller, interacting with gamma the same way as mu does for
    % the contractive PnP controller, hence---
    omega=mu;
end

%% Simulation

TSPAN = [0 flowTime];
JSPAN = [0 jumpTime];

% Rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 2;

options = odeset('RelTol',1e-6,'MaxStep',maxStep); 

% Simulate
[t,j,x] = HyEQsolver(@FlowMap,@JumpMap,@FlowSet,@JumpSet,x0Alt,TSPAN,JSPAN,rule,options,'ode45');
%[t,x] = ode45(@FlowMap,TSPAN,x0Alt');
%[t,x] = ode23(@FlowMap,TSPAN,x0Alt');

%% Weed out the chaff: remove redundant data frames to speed things up
x=permute(reshape(x,[length(t),stateDim,agentNum]),[3,2,1]);

% Remove unnecessary data frames by sampling time at near-uniform intervals
redundantFrames=[];
frameNum=1;
for step = 2:length(t)
    if floor(t(step-1)/timeBtwFrames)<floor(t(step)/timeBtwFrames) && frameNum<Nframes
        frameNum=frameNum+1;
    else
        redundantFrames=[redundantFrames,step];
    end
end
t(redundantFrames)=[];
j(redundantFrames)=[];
x(:,:,redundantFrames)=[];

save('simData.mat');

%% Prepare Main Figure

figHeight=800;
topViewWidth=700;
edgeLengthsWidth=700;

video_fig=figure('Position',[0,0,topViewWidth+edgeLengthsWidth,figHeight]);
%pnp_ax=subplot(1,2,1);
pnp_ax=subplot('Position',[0.05,0.1,0.5,0.8]);
%el_ax=subplot(1,2,2);
el_ax=subplot('Position',[0.6,0.1,0.35,0.8]);

% Size the topView plot
daspect(pnp_ax,[1,1,1]);

axes(pnp_ax);
axis equal
hold(pnp_ax,'on');

% Plot workspace
obsPlot=plot(workspacePolygon,'FaceColor','black','FaceAlpha',0.05,'DisplayName','$\quad$ Workspace ($\Omega$)');
if ~workspaceIsSphereWorld
    obsBufferPlot=plot(obstacleBuffer,'EdgeColor','none','FaceColor','red','FaceAlpha',0.1,'DisplayName','$\quad$ Obstacle buffers');
end
obsCirclesPlot=plot(obsDisks,'LineStyle','--','EdgeColor','red','FaceColor','red','FaceAlpha',0,'DisplayName','$\quad$ Obstacle internal radii');
% Plot the target, obstacle centers, obstacle internal radii
targetPlot=plot(pnp_ax,xStar(1),xStar(2),'*','DisplayName',"$\quad$ Target ($x^{\ast}$)");
centersPlot=plot(pnp_ax,obstacleCenters(:,1),obstacleCenters(:,2),'*','DisplayName','$\quad$ Obstacle Centers ($x_j^{\ast}$)');
% Plot the graph
masGraph=commGraph(x0); %prepare the graph using initial condition
nodeLabels=zeros(agentNum,1); %generate node labels, marking the leader
for j=1:agentNum
    nodeLabels(j,1)=j;
end
nodeLabels=string(nodeLabels);
nodeLabels(leader)=string(leader)+'L';
masGraphPlot=plot(pnp_ax,masGraph,'XData',x0(:,1),'YData',x0(:,2),'DisplayName',"$\quad$ Communication graph ($\mathcal{G}$)",'EdgeColor','#A2142F','Nodecolor','#A2142F','Marker','o','Markersize',6,'NodeLabel',nodeLabels); %generate the graph plot

% Plot navigation field to the target
xg=wXmin+0.5:1:wXmax-0.5; yg=wYmin+0.5:1:wYmax-0.5;
[X,Y]=meshgrid(xg,yg);
L=length(X);
U=zeros(L,L);
V=zeros(L,L);
for idx=1:L
    for idy=1:L
        v=navfq([X(idx,idy);Y(idx,idy)],workspacePolygon);
        U(idx,idy)=v(1,1);
        V(idx,idy)=v(2,1);
    end
end
navfield=quiver(pnp_ax,X,Y,U,V,'DisplayName','$\quad$ Navigation Field to $x^\ast$','Color','#4DBEEE');

% Prepare title,legend,etc...
if softSpring==0
    pnpGainString=", $\mu=$"+string(mu);
    pnpTypeString="contractive";
else
    pnpGainString=", $\omega=$"+string(mu);
    pnpTypeString="lazy";
end
if independentNav==1
    titleString="Agents navigate independently of each other";
    basic_info="$t=$";
else
    titleString="MAS motion under the "+pnpTypeString+" PnP controller";
    basic_info="R="+string(R)+pnpGainString+", $\gamma=$"+string(gamma)+", $t=$";
end
title(pnp_ax,titleString,'Fontsize',24);
subtitle(pnp_ax,basic_info+string(round(0.,3)),'Interpreter','latex','FontSize',16);
legend(pnp_ax,'Interpreter','latex','FontSize',16,'Location','northwest');

hold(pnp_ax,'off')

%% Generate edgeLengths subfigure
edgeLength = zeros(length(t),agentNum-1);
for i = 1:agentNum
    for step = 1:length(t)
        if i==agentNum
            edgeLength(step,i) = norm(x(i,:,step)-x(1,:,step))/R;
        else
            edgeLength(step,i) = norm(x(i,:,step)-x(i+1,:,step))/R;
        end
    end
end

hold(el_ax,'on')
plot(el_ax,t,varrho/R*ones(1,length(t)),'LineStyle','--','LineWidth',1.5,'Color','blue','DisplayName','$d=\varrho/R\qquad$');
plot(el_ax,t,varrhostar/R*ones(1,length(t)),'LineStyle','--','LineWidth',1.5,'Color','red','DisplayName','$d=\varrho^\ast/R\qquad$');
for i = 1:agentNum
  plot(el_ax,t,edgeLength(:,i),'LineWidth',0.5*(1+agentNum-i),'DisplayName','$i=$'+string(i));
end
eL_plot=plot(el_ax,t(1)*ones(1,agentNum),edgeLength(1,:),'LineWidth',0.5,'Color','b','Marker','o','DisplayName','$t=$'+string(round(t(1),3)));
legend(el_ax,'Interpreter','latex','FontSize',16,'Location','southeast');
if independentNav==1
    ylim(el_ax,[0,5]);
else
    ylim(el_ax,[0,1.2]);
end
xlabel(el_ax,'$t$\,=\,Time elapsed','Interpreter','latex','FontSize',16);
ylabel(el_ax,'$d=\Vert x_i - x_{i+1} \Vert/R$','Interpreter','latex','FontSize',16);

hold(el_ax,'off')

%% Generate the separate edge lengths figure

el_fig=figure('Position',[100,100,100+edgeLengthsWidth,50+figHeight]);
el_ax2=copyobj(el_ax,el_fig);
el_ax2.Position=[0.1,0.1,0.875,0.85];
hold(el_ax2,'on')
el_ax2.XAxis.FontSize=16;
el_ax2.YAxis.FontSize=16;
xlabel(el_ax2,'$t$\,=\,Time elapsed','Interpreter','latex','FontSize',20);
ylabel(el_ax2,'$d=\Vert x_i - x_{i+1} \Vert/R$','Interpreter','latex','FontSize',20);
legend(el_ax2,'Interpreter','latex','FontSize',16,'Location','southeast');
hold(el_ax2,'off')
imwrite(frame2im(getframe(el_fig)),strcat("edge-lengths-",string(1),".jpg"),'jpeg');

%% Generate the tracks figure

tracks_fig=figure('Position',[200,200,100+topViewWidth,50+figHeight]);

tracks_ax=copyobj(pnp_ax,tracks_fig);
tracks_ax.Position=[0.075,0.05,0.875,0.85];
hold(tracks_ax,'on')
tracks_ax.XAxis.FontSize=16;
tracks_ax.YAxis.FontSize=16;
title(tracks_ax,titleString,'Fontsize',24);
subtitle(tracks_ax,basic_info+string(round(0.,3)),'Interpreter','latex','FontSize',20);
legend(tracks_ax,'Interpreter','latex','FontSize',16,'Location','northwest');

% draw individual agent tracks
for i=1:agentNum
    xCoords=squeeze(x(i,1,:));
    yCoords=squeeze(x(i,2,:));
    track=plot(tracks_ax,xCoords,yCoords,'LineWidth',2);
    track.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% release and display image
hold(tracks_ax,'off')

% save the image
imwrite(frame2im(getframe(tracks_fig)),strcat("pnp_tracks-",string(1),".jpg"),'jpeg');


%% Generate the fronts figure

fronts_fig=figure('Position',[300,300,topViewWidth,figHeight]);

fronts_ax=copyobj(pnp_ax,fronts_fig);
fronts_ax.Position=[0.05,0.1,0.9,0.8];
hold(fronts_ax,'on')

legend(fronts_ax,'Interpreter','latex','FontSize',16,'Location','southeast');

% draw evolving MAS fronts
for step=100:100:length(t)
    xCoords=squeeze(x(:,1,step));
    yCoords=squeeze(x(:,2,step));
    pnpFront=plot(fronts_ax,xCoords,yCoords,'LineWidth',3);
    pnpFront.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% release and display image
hold(fronts_ax,'off')

% save the image
imwrite(frame2im(getframe(fronts_fig)),strcat("pnp_fronts-",string(1),".jpg"),'jpeg');


%% Generate evolving graphs video

% video parameters
pnp_movie(length(t))=struct('cdata',[],'colormap',[]);

% calculate frame (viewport) to be captured
rect=video_fig.Position;

%Capture the first frame and save it as a JPG
%video_fig.Visible='off';
drawnow
pnp_movie(1)=getframe(video_fig);%,rect);
imwrite(frame2im(pnp_movie(1)),strcat("pnp_movie_poster-",string(1),".jpg"),'jpeg');
%Capture the rest of the movie
for frame = 2:length(t)
    if mod(frame,10)==0
        frame %output the frame number to report progress
    end
    %tG=commGraph(x(:,:,frame);
    subtitle(pnp_ax,basic_info+string(round(t(frame),3)),'Interpreter','latex');
    masGraphPlot.XData=x(:,1,frame);
    masGraphPlot.YData=x(:,2,frame);
    eL_plot.XData=t(frame)*ones(1,agentNum);
    eL_plot.YData=edgeLength(frame,:);
    eL_plot.DisplayName='$t=$'+string(round(t(frame),3));
    drawnow
    pnp_movie(frame)=getframe(video_fig);%,rect);
    if WantSplitFrames==1
            imwrite(frame2im(pnp_movie(frame)),strcat("./split/pnp_movie_poster-",string(frame),".jpg"),'jpeg');
    end
end
%video_fig.Visible='on';

%Save the movie as MP4
if WantMovie==1
    %movie_fig=figure('Position',[0,0,topViewWidth+edgeLengthsWidth,figHeight]);
    %movie(movie_fig,pnp_movie,1,frameRate)
    myWriter=VideoWriter('pnp_movie','MPEG-4');
    myWriter.FrameRate=frameRate;   
    myWriter.Quality=100;
    open(myWriter);
    writeVideo(myWriter,pnp_movie);
    close(myWriter);
end

%% Function and Set Definitions

function cG = commGraph(conf)
    global agentNum R adjMat
    dM = zeros(agentNum,agentNum);
    for i = 1:agentNum
        for j = 1:agentNum
            if (norm(conf(i,:)-conf(j,:)) <= R) && (i ~= j)
                dM(i,j) = adjMat(i,j);
            end
        end
    end
    cG = graph(dM);
end

function v  = FlowSet(x) 
    v = 1;
end

function v  = JumpSet(x) 
    v = 0;
end


function x_dot = FlowMap(x,t)
    global independentNav agentNum stateDim adjMat gamma xStar leader
    x_dot = zeros(size(x));
    if independentNav==1
        % if every agent navigates independently:
        for i = 1:agentNum
           index_i = (1:stateDim) + stateDim*ones(1,stateDim)*(i-1);
           x_dot(index_i)=gamma*frakN(xStar,x(index_i));
        end
    else
        % if the PnP cooperative controller is applied, compute it:
        for i = 1:agentNum
            %compute the contribution of the PnP field for agent i:
            index_i = (1:stateDim) + stateDim*ones(1,stateDim)*(i-1);
            temp = zeros(stateDim,1);
            for j = 1:agentNum
                if adjMat(i,j) == 1
                    index_j = (1:stateDim) + stateDim*ones(1,stateDim)*(j-1);
                    temp = temp + edgeWeight(x(index_j),x(index_i))*frakN(x(index_j),x(index_i)); 
                end
                if j == agentNum
                    x_dot(index_i) = temp;
                end 
            end
        end
    
        %compute the task component V:
        V = zeros(agentNum*stateDim,1); 
        i = leader; % The leader is agent 1 
        index_i = (1:stateDim) + stateDim*ones(1,stateDim)*(i-1);
        V(index_i) = gamma*frakN(xStar,x(index_i))-x_dot(index_i); % Task component for asymmetric leader      
        %V(index_i) = gamma*frakN(xStar,x(index_i)); % Task component for symmetric leader      
    
        %augment the PnP field with the task component:
        x_dot = x_dot + V;
    end
end

function x_plus = JumpMap(x)
    x_plus = x;
end

function boo=adjQ(y,z)
    global R
    if norm(y-z)<=R
        boo=1;
    else
        boo=0;
    end
end


%% edge tension, edge weights

function output = edgeWeight(y,z)
    % Computes the assymmetric edge weights of the PnP controller
    s = norm(y-z);
    output = (r(s)*s^2)/(frakN(y,z)'*(y-z));
end

function output = r(s)
    % Edge tension function
    global softSpring R varrho mu omega alpha
    if s <= varrho
        output = (1-softSpring)*mu;
    elseif s <= R && s>= varrho
        output = (1-softSpring)*mu+omega*(s-varrho)^(1+alpha);
    else
        output = 0; %zeroing out results in the controller "breaking" edges of length greater than R.
    end
end

%% Navigation field

function v=navfq(pt,poly)
    % returns the navigation field to the target if pt is inside wsp;
    % otherwise returns a NaN vector
    % poly is a polyshape object
    % pt is a stateDim*1 vector
    global xStar stateDim
    if isinterior(poly,pt.')
        v=frakN(xStar,pt);
    else
        v=NaN(stateDim,1,'double');
    end
end

% Update this function with new options (obstacle types etc.)
function n=frakN(y,z)
    % Returns the normalized navigation field, depending on the user's
    % choice of environment / type of field
    % 
    global workspaceIsSphereWorld stateDim
    if workspaceIsSphereWorld
        n=frakN_sphere(y,z);
    else
        n=frakN_star(y,z);
    end
    if norm(n)<10^(-6)
        n=zeros(stateDim,1);
    else
        n=n/norm(n);
    end

end

%% Navigation field: Sphere world

function n=frakN_sphere(y,z)
    % Computes the navigation field by solving the associated quadratic
    % program
    global stateDim workspaceMatrix workspaceCoefficients
    A=[workspaceMatrix; safetyMatrix(z)];
    b=[workspaceCoefficients-workspaceMatrix*y; safetyCoefficients(y,z)];
    options=optimoptions('quadprog','Algorithm','active-set','Display','off');
    n=y+quadprog(eye(stateDim),zeros(stateDim,1),A,b,[],[],[],[],z-y,options)-z;
end

function m=safetyMatrix(z)
    % Computes the coefficient matrix describing the safe polytope at the point z
    % -  z is entered as a column vector
    global stateDim obstacleNum obstacleCenters
    m=zeros(obstacleNum,stateDim);
    for i=1:obstacleNum
        % compute the i-th row of the constraint matrix
        m(i,:)=obstacleCenters(i,:)-z';
    end
end

function c=obsDists(z)
    % Computes the column vector of distances of z to the obstacle centers
    % -  z is entered as a column vector
    global obstacleNum obstacleCenters
    c=zeros(obstacleNum,1);
    for i=1:obstacleNum
        col=z-obstacleCenters(i,:)';
        c(i,1)=c(i,1)+sqrt((col')*col);
    end
end

function b=safetyCoefficients(y,z)
    % Computes the column vector of inhomogeneous coefficients for the
    % translated quadratic optimization problem.
    % -  y,z are entered as a column vectors
    global obstacleNum obstacleRadii
    b=zeros(obstacleNum,1);
    dists=obsDists(z);
    cons=safetyMatrix(z);
    b=b+0.5*(dists.*dists-obstacleRadii.*dists)+cons*(z-y);
end

%% Navigation field: distorted sphere world
function n=frakN_star(y,z)
    global stateDim
    newy=workspaceToSphereWorld(y);
    newz=workspaceToSphereWorld(z);
    A=det(workspaceToSphereWorldDeriv(z));
    n=inv(workspaceToSphereWorldDeriv(z))*frakN_sphere(newy,newz);
end

function val=barrier(j,p)
    % Define the star-shaped obstacles. Must make sure that:
    % 1. The j-th obstacle is star-shaped relative to obstacleCenters(j,1);
    % 2. The j-th obstacle contains the ball of radius obstacleRadii(j,1)
    %    about the point obstacleCenters(j,1);
    % 3. The sets [barrier(j,p)<=obstacleClearances(j,1)] are pairwise
    %    disjoint.
    global stateDim obstacleNum obstacleCenters obstacleRadii obstacleClearance
    oc=obstacleCenters(j,:)';
    obsrad=obstacleRadii(j);
    relpos=p-oc;
    [theta,rho]=cart2pol(relpos(1),relpos(2));
    r0=barrierCurve(j,theta);
    val=rho^2-r0^2;
    %val=obsrad*(rho/r0-1);
end

function vec=barrierDeriv(j,p)
    % compute the derivative (=transpose of the gradient) of the jth
    % barrier at the point p
    global stateDim obstacleNum obstacleCenters obstacleRadii obstacleClearance
    oc=obstacleCenters(j,:);
    obsrad=obstacleRadii(j);
    relpos=p'-oc;
    [theta,rho]=cart2pol(relpos(1),relpos(2));
    u=[cos(theta),sin(theta)];
    v=[sin(theta),-cos(theta)];
    r0=barrierCurve(j,theta);
    r1=barrierCurveDeriv(j,theta);
    vec=2*(relpos+(r0*r1/rho)*v);
    %vec=-obsrad/(r0^2)*(r0*u+r1*v);
end

function out=barrierCurve(j,th)
    global obstacleRadii
    %syms th radius(th)
    obsrad=obstacleRadii(j);
    switch j
        case 2 %first obstacle boundary
            Npeaks=3;
            out=1.1*obsrad+(1+sin(Npeaks*th));
        case 1 %second obstacle boundary
            Npeaks=7;
            out=1.1*obsrad+(1+sin(Npeaks*th));
        case 3 %third obstacle boundary
            offsetAngle=pi/6;
            Npeaks=2;
            out=1.1*obsrad+2*(1+cos(Npeaks*(th+offsetAngle)));
    end
end

function out=barrierCurveDeriv(j,th)
    global obstacleRadii
    %syms th radius(th)
    obsrad=obstacleRadii(j);
    switch j
        case 2 %first obstacle boundary
            Npeaks=3;
            out=Npeaks*cos(Npeaks*th);
        case 1 %second obstacle boundary
            Npeaks=7;
            out=Npeaks*cos(Npeaks*th);
        case 3 %third obstacle boundary
            offsetAngle=pi/6;
            Npeaks=2;
            out=-2*Npeaks*sin(Npeaks*(th+offsetAngle));
    end
end


function y=stitch(j,x)
    % stitching function
    global obstacleClearance
    epsilon=obstacleClearance(j,1);
    y=bump(epsilon-x)/bump(epsilon);
end

function y=stitchDeriv(j,x)
    % derivative of the stitching function
    global obstacleClearance
    epsilon=obstacleClearance(j,1);
    y=-bumpDeriv(epsilon-x)/bump(epsilon);
end

function v=stitchingMap(p)
    % compute the vector of stitching values $\sigma_i$ at p, including the
    % neutral value, $\sigma_d$, provided in the last coordinate of the
    % output.
    global obstacleNum
    v=zeros(obstacleNum+1,1);
    v(obstacleNum+1,1)=1;
    for j=1:obstacleNum
        v(j,1)=v(j,1)+stitch(j,barrier(j,p));
        v(obstacleNum+1,1)=v(obstacleNum+1,1)-v(j,1);
    end
end

function smd=stitchingMapDeriv(p)
    % compute the (obstacleNum+1)*stateDim derivative of the stitching map
    global stateDim obstacleNum
    smd=zeros(obstacleNum+1,stateDim);
    for j=1:obstacleNum
        smd(j,:)=(stitchDeriv(j,barrier(j,p)))*barrierDeriv(j,p);
    end
    smd(j+1,:)=-sum(smd);
end

function q=workspaceToSphereWorld(p)
    % compute the transformation of the workspace to the associated sphere
    % world. This ASSUMES a barrier function barrier(j,-) is provided for
    % each obstacle j=1,...,obstacleNum; that each obstacle O_j is star-shaped
    % relative to the point obstacleCenters(j,1); and that each obstacle
    % O_j contains the ball about obstaclesCenters(j,1) with radius
    % obstacleRadii(j,1).
    global obstacleNum obstacleCenters obstacleRadii obstacleClearance
    stitches=stitchingMap(p);
    q=stitches(obstacleNum+1,1)*p;
    for j=1:obstacleNum
        oc=obstacleCenters(j,:)';
        rad=obstacleRadii(j,1);
        stitch=stitches(j,1);
        relpos=p-oc;
        relnorm=norm(relpos);
        q=q+stitch*(oc+rad*relpos/relnorm);
    end
end

function D=workspaceToSphereWorldDeriv(p)
    % compute the derivative of the transformation of the workspace to the 
    % associated sphere world (see workspaceToSphereWorld above).
    global stateDim obstacleNum obstacleCenters obstacleRadii
    stitches=stitchingMap(p);
    stitchesDeriv=stitchingMapDeriv(p);
    D=stitches(obstacleNum+1)*eye(stateDim);
    for j=1:obstacleNum
        oc=obstacleCenters(j,:)';
        rad=obstacleRadii(j,1);
        stitch=stitches(j,1);
        relpos=p-oc;
        relnorm=norm(relpos);
        D=D+(stitch*rad/relnorm)*eye(stateDim);
        D=D+((rad/relnorm)-1)*relpos*stitchesDeriv(j,:);
        D=D-stitch*rad*(relnorm^(-3))*(relpos*relpos');
    end
end

%% Bump functions
function y=bump(x)
    y=cInftyBump(x);
end

function y=bumpDeriv(x)
    y=cInftyBumpDeriv(x);
end

function y=cInftyBump(x)
    % C^\infty bump function
    if x<=0
        y=0;
    else
        y=exp(-1./x);
    end
end

function y=cInftyBumpDeriv(x)
    %derivative of the C^\infty bump function
    if x<=0
        y=0;
    else
        y=exp(-1./x)*(x^(-2));
    end
end

function y=cOneBump(x)
    % C^1 bump function
    if x<=0
        y=0;
    else
        y=x^2;
    end
end

function y=cOneBumpDeriv(x)
    %Derivative of C^1 bump function
    if x<=0
        y=0;
    else
        y=2*x;
    end
end

function y=cTwoBump(x)
    % C^1 bump function
    if x<=0
        y=0;
    else
        y=x^3;
    end
end

function y=cTwoBumpDeriv(x)
    %Derivative of C^1 bump function
    if x<=0
        y=0;
    else
        y=3*x^2;
    end
end


