%MATLAB CODE
%EVALUATES THE STATE DISTRIBUTION FOR THE
%EXACT SYSTEM USING THE POWER METHOD.

%AUTHOR: ANDERS REENBERG ANDERSEN

%------------------------
%   PREAMPLE
%------------------------

clear all; clc;

%number of CPU cores
parWorkers = 4;

%------------------------
%   PARAMETERS
%------------------------

%label
lab_test = "SomeLabel";

%rho
rutil = 0.80;

%number of wards
nWards = 3;
%capacity
capacity = [3,3,3];
%service rates
serRate_all = 0.1;
mu = [serRate_all,serRate_all,serRate_all];
%arrival rates
arrRate_all = rutil*(serRate_all*capacity(1));
lambda = [arrRate_all,arrRate_all,arrRate_all];


%routing matrix
R = [0.0 0.5 0.5;
     0.5 0.0 0.5;
     0.5 0.5 0.0];
 
%state space size
n = getStateSpaceSize(capacity,nWards);
%generate state space (for small models)
S = getStateSpace(capacity,nWards);
%transition rate matrix (for small models)
Q = zeros(n,n); 

%------------------------
%   RUN
%------------------------

%calculate off-diagonal elements
parfor (i = 1:n,parWorkers) 
    for j = 1:n
        if i~=j
            Q(i,j) = getTransitionRate(i,j,mu,lambda,capacity,S,R,nWards);
        end
    end
end    

%calculate diagonal elements
qdiag = -sum(Q,2);
for i = 1:n
    Q(i,i) = qdiag(i); 
end

%save some variables
%save('Qmatrix','Q');
%csvwrite('Qmatrix.csv',Q);
%csvwrite('Svector.csv',S);

offset = 1e-3;
deltat = 1/max(abs(qdiag)) - offset;
%deltat = 1/max(abs(qdiag));
P = Q*deltat + eye(n);
Pt = transpose(P);

pi = rand(n,1);
sm = sum(pi);
pi_new = pi./sm;

m = 100;
tol = 1e-9;
eps = 1;
while eps>tol
    pi_old = pi_new;
    
    pi_new = Pt*pi_old;
    for i = 1:m
        pi_new = Pt*pi_new;
        %normalize
        sm = sum(pi_new);
        pi_new = pi_new./sm;
    end
        
    %evaluate convergence
    eps = max( (abs(pi_new-pi_old)./pi_new) );
    disp(eps)
end
pi = pi_new;

%marginal distributions
%for widx = 1:nWards
%    dist = marginalDist(pi,widx,nWards,capacity,S,n);
%    disp(dist)
%end

dist1 = marginalDist(pi,1,nWards,capacity,S,n);
dist2 = marginalDist(pi,2,nWards,capacity,S,n);
%dist3 = marginalDist(pi,3,nWards,capacity,S,n);

csvwrite(strcat(lab_test,"_Dist.csv"),[dist1 dist2])

%------------------------
%   FUNCTIONS
%------------------------

function dist = marginalDist(pi,widx,nWards,capacity,S,n)
    dist = zeros((capacity(widx)+1),1);
    for i = 1:n
        k = capUse(widx,S(i,:),nWards);
        idx = k+1;
        dist(idx) = dist(idx) + pi(i); 
    end    
end

function col = getColumn(widx,pdix,nWards)
    col = nWards*(widx-1) + pidx;
end

function [widx,pidx] = getIndices(col,nWards)
    widx = ceil(col/nWards);
    pidx = col-nWards*(widx-1);
end

function widx = getWard(col,nWards)
    widx = ceil(col/nWards);
end

function k = capUse(widx,s,nWards)
    startIdx = nWards*(widx-1)+1;
    stopIdx = startIdx+(nWards-1); 
    k = sum(s(startIdx:stopIdx));
end

function q = getTransitionRate(i,j,mu,lambda,capacity,S,R,nWards)
    diff = S(j,:)-S(i,:);
    ch = length(find(diff~=0)); %number of state changes
    
    q=0;
    if ch == 1
        arr = length(find(diff==1));
        dis = length(find(diff==-1));
    
        if arr==1 && dis==0 %arrival
            cidx = find(diff==1);
            [widx,pidx] = getIndices(cidx,nWards);
            if widx==pidx
                q = lambda(widx);
            elseif capUse(pidx,S(i,:),nWards)==capacity(pidx)
                q = lambda(pidx)*R(pidx,widx);
            end    
        elseif arr==0 && dis==1 %discharge    
            cidx = find(diff==-1);
            [widx,pidx] = getIndices(cidx,nWards);
            q = S(i,cidx)*mu(pidx);
        end
    end
end

function n = getStateSpaceSize(capacity,nWards)
    z = 1/factorial(nWards);
    n = 1;
    for j = 1:nWards
        w = z;
        for i = 1:nWards
            w = w*(capacity(j)+i);
        end
        n = n*w;
    end
end

function S = getStateSpace(capacity,nWards)
    n = getStateSpaceSize(capacity,nWards);
    sz = nWards^2; %state size
    S = zeros(n,sz);
    Kuse = zeros(1,nWards);
    
    for i = 2:n
        S(i,:) = S(i-1,:);
        j = sz;
        while j > -1
            widx = getWard(j,nWards);
            if Kuse(widx)<capacity(widx)
                S(i,j) = S(i,j) + 1;
                Kuse(widx) = Kuse(widx)+1;
                j = -1;
            elseif S(i,j)>0
                Kuse(widx) = Kuse(widx) - S(i,j);
                S(i,j) = 0;
                j = j - 1;
            else
                j = j - 1;
            end
        end
        
    end
    
end

