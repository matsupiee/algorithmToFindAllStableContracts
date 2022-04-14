%% setup 

clear all;

global M W N_m N_w N_t

% Parameters
N_w = 5; % # of women
N_m = 6; % # of men
N_t = 2; % (# of contracts)/2

% Utility matrix
% set seed
rng(6)
M = randn(N_m,N_w); % Men's utility over women
W = randn(N_w,N_m); % Women's utility over men

% Create contract space 
dif = 0.1;
TR = zeros(N_t,1); % matrix of contracts
min_t = 0;

for n=1:N_t
    % positive payments
    TR(n,1) = (exp(min_t + (n-1)*dif) + exp(min_t + n*dif))/2;
end

TR = [-flipud(TR);TR]; % contracts (allow both pos/neg payments); 2*N_t

TR_expand_m = repmat(TR',N_m,N_w); % N_w x (N_m x 2N_t)
TR_expand_m = sort(TR_expand_m,2,'ascend');
P_m = [zeros(N_m,N_w*N_t),repmat(M,1,N_t)] - 0.1*TR_expand_m; % men's utility over men with transfers
rng(100) % set seed for tiebreaking
tiebreak_m = 0.001 * randn(size(P_m)); % previously 0.001*[repmat(M,1,N_t),zeros(N_m,N_w*N_t)];
M = P_m + tiebreak_m;

TR_expand_w = repmat(TR',N_w,N_m);
TR_expand_w = sort(TR_expand_w,2,'ascend'); 
P_w = [repmat(W,1,N_t),zeros(N_w,N_m*N_t)] + 0.1*TR_expand_w; % women's utility over men with transfers
rng(10) % set seed for tiebreaking
tiebreak_w = 0.001 * randn(size(P_w)); % previously 0.001*[zeros(N_w,N_m*N_t),repmat(W,1,N_t)];
W = P_w + tiebreak_w;

M=M';
W=W';
% (i,j)th element is j's utility from i-th contract

%% prematching
% construct largest prematching (lv) for man (lvm) and woman (lvw)

% lvm
[U_lvm,lvm]=max(M);
% if no one is acceptable for i, matche i with empty set (0)
for i=1:length(lvm)
    if U_lvm(1,i)<0
        lvm(1,i)=0;
    end
end

% lvw
[U_lvw,lvw]=max(W);
for i=1:length(lvw)
    if U_lvw(1,i)<0
        lvw(1,i)=0;
    end
end

% largest prematching (lv)
lv=[lvm,lvw];

%% iteration

% T2 iteration on largest prematching v (lv)
T2lv_old = lv;
T2lv_new = functionT_with_contract(M,W,functionT_with_contract(M,W,T2lv_old));

while isequal(T2lv_new,T2lv_old) == 0
    T2lv_old = T2lv_new;
    T2lv_new = functionT_with_contract(M,W,functionT_with_contract(M,W,T2lv_old));
end

% largest fixed-point prematching
lv_T2fix = T2lv_new;

% T^2 iteration on smallest prematching v (v0, since all elements are zero)
T2v0_old = zeros(1,N_m+N_w);
T2v0_new = functionT_with_contract(M,W,functionT_with_contract(M,W,T2v0_old));

while isequal(T2v0_new,T2v0_old) == 0
    T2v0_old = T2v0_new;
    T2v0_new = functionT_with_contract(M,W,functionT_with_contract(M,W,T2v0_old));
end

% smallest fixed-point prematching
v0_T2fix = T2v0_new;

% algorithm
if isequal(lv_T2fix,functionT_with_contract(M,W,lv_T2fix))==1
    allcore=lv_T2fix;
    
else 
    % subroutine
    
    % apply subroutine to Q={lv}
    [allcore,allQ]=subroutine_step_with_contract(lv_T2fix,v0_T2fix);
    
    % repeat until Q becomes empty
    while isequal(allQ,zeros(0,N_m+N_w)) == 0
        
        % index for prematching in Q
        x=size(allQ,1);
        
        % storing Q in each loop
        store_Q=zeros(0,N_m+N_w);
        
        % do subroutine steps for each prematching in Q
        for k=1:x
            v=allQ(k,:);
            [core_k,Q_k]=subroutine_step_with_contract(v,v0_T2fix);
            allcore =[allcore;core_k];
            store_Q =[store_Q;Q_k];
        end
        
        % core and Q
        allcore = unique(allcore,'stable','rows');
        allQ = unique(store_Q,'stable','rows');
        
       
    end     
    
end


%% cores


% check cores
temp_allcore_num = height(allcore);
core_check = ones(1,temp_allcore_num);
for i = 1:temp_allcore_num
    x = stableOrNot_with_contract(allcore(i,:));
    if x ~= 1
        core_check(1,i) = x;
    end
end


% divide cores into men's part and women's part
allcore_m =allcore(:,1:N_m);
allcore_w =allcore(:,N_m+1:N_m+N_w);

% number of stable matchings
allcore_num = height(allcore); 

% matching mate in core
allcore_m_agent = allcore_m;
allcore_w_agent = allcore_w;

for i=1:N_m
    for j=1:allcore_num
        if allcore_m(j,i)==0
            allcore_m_agent(j,i)=0;
        elseif mod(allcore_m(j,i),N_w)==0
            allcore_m_agent(j,i)=N_w;
        else
            allcore_m_agent(j,i)=mod(allcore_m(j,i),N_w);
        end
    end
end
for i=1:N_w
    for j=1:allcore_num
        if allcore_w(j,i)==0
            allcore_w_agent(j,i)=0;
        elseif mod(allcore_w(j,i),N_m)==0
            allcore_w_agent(j,i)=N_m;
        else
            allcore_w_agent(j,i)=mod(allcore_w(j,i),N_m);
        end
    end
end
allcore_agent = [allcore_m_agent,allcore_w_agent];

