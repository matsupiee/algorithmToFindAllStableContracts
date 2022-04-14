function[core,Q]=subroutine_step_with_contract(v,sf)

global M W N_m N_w N_t

% divide prematching v into men's part and women's part
vm=v(1,1:N_m);
vw=v(1,N_m+1:end);

% divide the smallest prematching sf into men's part and women's part
sfm=sf(1,1:N_m);
sfw=sf(1,N_m+1:end);

%まず、各エージェントについて、prematching v以上に好ましい相手とマッチしたときの利得が0である
%選好modifiedMとmodifiedWを作る。
% Make modified (truncated) preferences mM and mW so that each agent has 0 payoff
% from the better mates than his/her mate in the prematching v
mM = M;
mW = W;


% modified preference of men mM
for i=1:N_m
    
    for j=1:N_w*2*N_t
        
        %もし、prematching v　において、man iのマッチ先が空集合であれば、それより大きいマッチ相手の利得も0。
        % If man i is matched to the empty set in prematching v, contracts whose utility are larger than 0 give 0 payoff 
        if vm(1,i)==0
            if mM(j,i) > 0
                mM(j,i)=0;
            end
            
        % If man i is matched by contract j in prematching v, all better contracts give 0 payoff 
        elseif M(vm(1,i),i) < M(j,i)
            mM(j,i)=0;
            
        end
    end
end



% modified preference of women mW
for i=1:N_w
    
    for j=1:N_m*2*N_t
        
        % If woman i is matched to the empty set in prematching v, all contracts give 0 payoff 
        if vw(1,i)==0
            if mW(j,i)>0
                mW(j,i)=0;
            end
        % If woman i is matched by contract j in prematching v, all better contracts give 0 payoff 
        elseif W(vw(1,i),i) < W(j,i)
            mW(j,i)=0;
        end
    end
end          


% matrices for the results
core=zeros(0,N_m+N_w); % matrix storing stable matchings
Q=zeros(0,N_m+N_w); % matrix storing prematching used in the next subroutine




%man i について、vm(1,i)がsfm(1,i)よりも大きければ、step1&2を回す。
%prematching vにおいてiのマッチ先が空集合でない必要がある。
% If vm(1,i) > sfm(1,i) for man i, proceed to step 1 & 2.
% This is equivalent to i(f,v)R(f)v_(f) under strong preference.
% man i should not be matched to the empty set in prematching v.

for i=1:N_m
    
    % if man i is matched to the empty set in sfm
    if sfm(1,i)== 0
        
        % vm(1,i)>0 to satisfy the condition
        if vm(1,i)>0
            
            % modified preference for man i
            mMi=mM;
            mMi(vm(1,i),i)=0;
            [C,vmi]=max(mMi);
            
            % if in column j the largest value <= 0,
            %j will choose outside option.
            for j=1:N_m
                if C(1,j)<=0
                    vmi(1,j)=0;
                end
            end
            
            % largest prematching vi
            vi=[vmi,vw];
            
            % T^2 iteration on vi
            T2vi_old = vi;
            T2vi_new = functionT_with_contract(mMi,mW,functionT_with_contract(mMi,mW,T2vi_old));
            while isequal(T2vi_new,T2vi_old) == 0
                 T2vi_old = T2vi_new;
                 T2vi_new = functionT_with_contract(mMi,mW,functionT_with_contract(mMi,mW,T2vi_old));
            end
            
            % largest fixed-point prematching
            vi_T2fix = T2vi_new;
                        
            % If Tvi_T2fix = vi_T2fix, add vi_T2fix to the core
            Tvi_T2fix = functionT_with_contract(mMi,mW,vi_T2fix);
            
            if isequal(Tvi_T2fix,vi_T2fix)==1
                core=[core;vi_T2fix];
                
           
            % If vi_T2fix > sf,add vi_T2fix to Q
            elseif largerOrNot(mMi,mW,vi_T2fix,sf)==1
                Q=[Q;vi_T2fix];
                
            end
        end   
        
    % if man i is matched to some woman in sfm
    elseif vm(1,i)>0
        
        % If vm(1,i) > sfm(1,i) 
        if M(vm(1,i),i) > M(sfm(1,i),i)
            
            %iがprematching vにおけるマッチ先とマッチしたときの利得を0にしておく。
            % Let the payoff of i matched to the same woman as in the
            % prematching v be 0
            mMi=mM;
            mMi(vm(1,i),i)=0;
            [C,vmi]=max(mMi);
            
            % if in column j the largest value <= 0,
            % j will choose outside option.
            for j=1:N_m
                if C(1,j)<=0
                    vmi(1,j)=0;
                end
            end
            
            % largest prematching vi
            vi=[vmi,vw];
            
            % T^2 iteration on vi
            T2vi_old = vi;
            T2vi_new = functionT_with_contract(mMi,mW,functionT_with_contract(mMi,mW,T2vi_old));
            while isequal(T2vi_new,T2vi_old) == 0
                 T2vi_old = T2vi_new;
                 T2vi_new = functionT_with_contract(mMi,mW,functionT_with_contract(mMi,mW,T2vi_old));
            end
            
            % largest fixed-point prematching
            vi_T2fix = T2vi_new;
                        
            % If Tvi_T2fix = vi_T2fix, add vi_T2fix to the core
            Tvi_T2fix=functionT_with_contract(mMi,mW,vi_T2fix);
            
            if isequal(Tvi_T2fix,vi_T2fix)==1
                core=[core;vi_T2fix];
                
            % If vi_T2fix does not coincide with sf, vi_T2fix > sf
            elseif largerOrNot(mMi,mW,vi_T2fix,sf)==1
                Q=[Q;vi_T2fix];
                
            end
        end
    end
end



%woman i について、vw(1,i)がsfw(1,i)よりも大きければ、step1&2を回す。
% If vm(1,i) > sfm(1,i) for woman i, proceed to step 1 & 2.
% This is equivalent to i(f,v)R(f)v_(f) under strong preference.
% Woman i should not be matched to the empty set.

for i=1:N_w
    
    % if woman i is matched to the empty set in sfm
    if sfw(1,i)==0
        
        % vm(1,i)>0 to satisfy the condition
        if vw(1,i)>0
            
            % modified preference for woman i
            mWi=mW;
            mWi(vw(1,i),i)=0;
            [C,vwi]=max(mWi);
            
            % if in column j the largest value <= 0,
            % j will choose outside option.
            
            for j=1:N_w
                if C(1,j)<=0
                    vwi(1,j)=0;
                end
            end
            
            % largest prematching vi
            vi=[vm,vwi];
            
            % T^2 iteration on vi
            T2vi_old = vi;
            T2vi_new = functionT_with_contract(mM,mWi,functionT_with_contract(mM,mWi,T2vi_old));
            while isequal(T2vi_new,T2vi_old) == 0
                 T2vi_old = T2vi_new;
                 T2vi_new = functionT_with_contract(mM,mWi,functionT_with_contract(mM,mWi,T2vi_old));
            end
            
            % largest fixed-point prematching
            vi_T2fix = T2vi_new;
                        
            % If Tvi_T2fix = vi_T2fix, add vi_T2fix to the core
            Tvi_T2fix=functionT_with_contract(mM,mWi,vi_T2fix);
            
            if isequal(Tvi_T2fix,vi_T2fix)==1
                core=[core;vi_T2fix];
                
            % If vi_T2fix > sf, add it to Q
            elseif largerOrNot(mM,mWi,vi_T2fix,sf)==1
                Q=[Q;vi_T2fix];
                
            end
        end
        
    % if woman i is matched to some man in v
    elseif vw(1,i)>0
        
        % If vw(1,i) > sfw(1,i) 
        if W(vw(1,i),i) > W(sfw(1,i),i)
            
            % modified preference for woman i
            mWi=mW;
            mWi(vw(1,i),i)=0;
            [C,vwi]=max(mWi);
            
            % if in column j the largest value <= 0,
            % j will choose outside option.
            
            for j=1:N_w
                if C(1,j)<=0
                    vwi(1,j)=0;
                end
            end
            
            
            
            % largest prematching vi
            vi=[vm,vwi];
            
            % T^2 iteration on vi
            T2vi_old = vi;
            T2vi_new = functionT_with_contract(mM,mWi,functionT_with_contract(mM,mWi,T2vi_old));
            while isequal(T2vi_new,T2vi_old) == 0
                 T2vi_old = T2vi_new;
                 T2vi_new = functionT_with_contract(mM,mWi,functionT_with_contract(mM,mWi,T2vi_old));
            end
            
            % largest fixed-point prematching
            vi_T2fix = T2vi_new;
                        
            % If Tvi_T2fix = vi_T2fix, add vi_T2fix to the core
            Tvi_T2fix=functionT_with_contract(mM,mWi,vi_T2fix);
            
            if isequal(Tvi_T2fix,vi_T2fix)==1
                core=[core;vi_T2fix];
                
            % If vi_T2fix does not coincide with sf, vi_T2fix > sf
            elseif largerOrNot(mM,mWi,vi_T2fix,sf)==1
                Q=[Q;vi_T2fix];
                
            end
        end
    end
end
end





        
