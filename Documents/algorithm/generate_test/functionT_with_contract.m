function[Tv]=functionT_with_contract(M,W,v)

global N_m N_w N_t

% divide prematching into men's part and women's part
vm=v(1,1:N_m);
vw=v(1,N_m+1:end);

%Tv(man)について考える
%prematching vよりも好ましいコントラクトに1をつける行列を作る。
% Consider Tv(man)
% Matrix where contracts which each woman prefers to the contract in
% prematching v have 1
U=zeros(N_w*2*N_t,N_m);

for i=1:N_m
    
     for j = 1:N_w
         
         % if woman j is not matched to anyone (empty set)
         if vw(1,j) == 0
             
             % give 1 to a contract which gives positive utility of woman j
             for k = 1:2*N_t
                 x = k-1;
                 if W(N_m*x+i,j)>0
                     y = N_w*x+j;
                     U(y,i)=1;
                 end
             end
             
         % if woman j is matched to a man
         else
             
             % give 1 to a contract which gives larger utility of woman
             % than j's contract in the premaching
             for k = 1:2*N_t
                 x = k-1;
                 if W(N_m*x+i,j)>=W(vw(1,j),j)
                     y = N_w*x+j;
                     U(y,i)=1;
                 end
             end
         end
     end
end



%%Tv(woman)について考える。
% Consider Tv(woman)
% Matrix where contracts which each man prefers to the contract in
% prematching v have 1
V=zeros(N_m*2*N_t,N_w);

for i=1:N_w
    
     for j = 1:N_m
         
         % if man j is not matched to anyone (empty set)
         if vm(1,j) == 0
             
             for k = 1:2*N_t
                 x = k-1;
                 if M(N_w*x+i,j)>0
                     y = N_m*x+j;
                     V(y,i)=1;
                 end
             end
             
         % if man j is matched to a woman
         else
             
             for k = 1:2*N_t
                 x = k-1;
                 if M(N_w*x+i,j)>=M(vm(1,j),j)
                     y = N_m*x+j;
                     V(y,i)=1;
                 end
             end
         end
     end
end

% Each man/woman chooses the best contract from U/V
% For man i, put his utilities over the contracts which has 1 in U matrix.
%for i=1:N_m
%    for j=1:N_w*2*N_t
%        if U(j,i)==1
%            U(j,i)=M(j,i);
%        end
%    end   
% end

% (more efficient code)
U = U.*M;


%各列で最も利得が高い行のインデックスを返す
% index of the row which has the highest utility for each column
[Mx,vm]=max(U);

%全ての行が0の列では、1を返してしまうので、そのような列は0に変えておく。
% For an all-zero column, give index 0
for i=1:N_m
    if Mx(1,i)<=0
        vm(1,i)=0;
    end
end




% For woman i, put her utilities over the contracts which has 1 in V matrix.
%for i=1:N_w
%     for j = 1:N_m*2*N_t
%        if V(j,i)==1
%            V(j,i)=W(j,i);
%        end
%     end
%end
V = V.*W;

% index of the row which has the highest utility for each column
[Wx,vw]=max(V);

% For an all-zero column, give index 0
for i=1:N_w
    if Wx(1,i)<=0
        vw(1,i)=0;
    end
end
 
% Tv
Tv = [vm,vw];
end

 
 