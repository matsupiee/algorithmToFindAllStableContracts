function[x]=stableOrNot_with_contract(v)

global M W N_m N_w N_t

% If v is stable matching, return 1
x=1;

%% divide prematchings into men's part and women's part
vm=v(1,1:N_m);
vw=v(1,N_m+1:end);



%% feasibility check

%if agents are not matched to properly, return -1.
%elseif contract number is not correct, return -1.5

%man-side
for i=1:N_m
    
    % check whether man-i's mate has man i as her mate
    if vm(1,i)~=0         
        % If vm(1,i)=0, ignore man i
        
        % let man i's mate be iw
        iw=mod(vm(1,i),N_w);
        
        if iw ~=0 % iw is 1 ~ N_w-1
            
            % let y be the contract number for iw
            y=vw(1,iw);
            
            % if iw is not matched to anyone, not feasible
            if y == 0
                x = -1;
                return
                
            % if iw is matched to a man, check whether he is man i    
            % mod(y,N_m) is iw's mate
            elseif mod(y,N_m)~=i  
                
                % If iw is matched to man N_m, mod(y,N_m) is 0
                if mod(y,N_m)==0 
                    
                    % if i is not N_m, infeasible
                    if i ~= N_m
                        x = -1;
                        return
                    end
                    
                % otherwise, infeasible
                else
                    x = -1;
                    return
                end
            end
            
        % iw=0 means iw is woman N_w    
        else   
            iw = N_w;
            
            % let y be the contract number for iw
            y=vw(1,iw);
            
            % if iw is not matched to anyone, not feasible
            if y == 0
                x = -1;
                return
                
            % if iw is matched to a man, check whether he is man i    
            % mod(y,N_m) is iw's mat
            elseif mod(y,N_m)~=i  
                
                % If iw is matched to man N_m, mod(y,N_m) is 0
                if mod(y,N_m)==0  
                    if i ~= N_m
                        x = -1;
                        return
                    end
                    
                % otherwise, infeasible    
                else
                    x = -1;
                    return
                end
            end
        end
        
        
        % after checking man i and woman iw are matched properly,
        % check the contract number
        
        if x == 1
            ik=ceil(vm(1,i)/N_w); % man i's contract number 
            iwk=ceil(y/N_m); %woman iw's contract number
            if ik ~= iwk
                x = -1.5;
                return
            end
        end
        
    end
end

%woman-side
for i=1:N_w
    
    % if x = -1, finish
    if x == -1
        return
        
    % if not
    % ignore the case where woman i is matched to the empty set
    elseif vw(1,i)~=0 
        
        % let woman i's mate be im
        im=mod(vw(1,i),N_m); 
        
        % if im is matched to a woman, check whether she is woman i    
        if im ~=0 % im is 1 ~ N_m-1
            
            % let y be the contract number for im
            y=vm(1,im);  
            if y == 0
                x = -1;
                return

            % mod(y,N_w) is im's mate
            elseif mod(y,N_w)~=i 
                
                % If im is matched to woman N_w, mod(y,N_w) = 0
                if mod(y,N_w)==0  
                    if i ~= N_w
                        x = -1;
                        return
                    end
                else
                    x = -1;
                    return
                end
            end
            
        %im=0 means that im=N_m
        else   
            im = N_m;
            % im's contract
            y=vm(1,im);
            if y == 0
                x = -1;
                return
                
            
            elseif mod(y,N_w)~=i  
                if mod(y,N_w)==0  % y is a contract with N_w
                    if i ~= N_w
                        x = -1;
                        return
                    end
                else
                    x = -1;
                    return
                end
            end
        end
        
        % after checking woman i and man im are matched properly,
        % check the contract number
        if x == 1
            ik=ceil(vw(1,i)/N_m); %woman i's contract number
            imk=ceil(y/N_w); %man im's contract number
            if ik ~= imk
                x = -1.5;
            end
        end
    end
end


%% check individual rationality
% If not satisfied, x = -2
if x == 1
    for i=1:N_m
        if vm(1,i)~=0
            if M(vm(1,i),i)<0
                x=-2;
                return
            end
        end
    end
end

if x == 1
    for i=1:N_w
        if vw(1,i)~=0
            if W(vw(1,i),i)<0
                x=-2;
                return
            end
        end
    end
end

%% check stability

% if there is a blocking pair, return x = 0

if x==1
    for i=1:N_m
        
        % if man i is not matched to anyone
        if vm(1,i)==0 
            for j=1:N_w
                for k=1:2*N_t   
                    
                    % if man i prefers k-th contract with woman j to
                    % outside option
                    if M((k-1)*N_w+j,i)>0  
                        
                        % if woman j prefers k-th contract with man i to
                        % outside option
                        if W((k-1)*N_m+i,j)>0  
                            
                            if vw(1,j)==0
                                x=0;
                                return
                                
                            % if woman j is better off from matching v
                            elseif W(vw(1,j),j)<W((k-1)*N_m+i,j)
                                x=0;
                                return
                            end
                        end
                    end
                end
            end
            
            
            
        % if man i is matched to someone
        else    
            for j=1:N_w
                for k=1:2*N_t  
                    
                    % if man i prefers k-th contract with woman j to matching v
                    if M((k-1)*N_w+j,i)>M(vm(1,i),i)     
                        
                        % if woman j prefers k-th contract with man i to
                        % the outside option
                        if W((k-1)*N_m+i,j)>0  
                            
                            if vw(1,j)==0
                                x=0;
                                return
                            
                            % if woman j is better off from matching v
                            elseif W(vw(1,j),j)<W((k-1)*N_m+i,j)
                                x=0;
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end  



if x==1
    for i=1:N_w
        
        % if woman i is not matched to anyone
        if vw(1,i)==0  
            for j=1:N_m
                for k=1:2*N_t   
                    
                    % if woman i prefers k-th contract with man j to outside option
                    if W((k-1)*N_m+j,i)>0      
                        
                        % if man j prefers k-th contract with woman i to outside option
                        if M((k-1)*N_w+i,j)>0
                            
                            % if man j chooses the outside option in v,
                            % blocking
                            if vm(1,j)==0  
                                x=0;
                                return
                            
                            % if man j is better off from matching v
                            elseif M(vm(1,j),j)<M((k-1)*N_w+i,j)
                                x=0;
                                return
                            end
                        end
                    end
                end
            end
            
            
        % if woman i is matched to someone
        else    
            for j=1:N_m
                for k=1:2*N_t   
                    
                    % if woman i prefers k-th contract with man j to matching v
                    if W((k-1)*N_m+j,i)>W(vw(1,i),i)     
                        
                        % if man j prefers k-th contract with woman i to the outside option
                        if M((k-1)*N_w+i,j)>0  
                            if vm(1,j)==0
                                x=0;
                                return
                            elseif M(vm(1,j),j)<M((k-1)*N_w+i,j)
                                x=0;
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end                                       
                                    


