function[x]=largerOrNot(M,W,v,sf)

global N_m N_w

%% divide prematchings into men's part and women's part
vm=v(1,1:N_m);
vw=v(1,N_m+1:end);

sfm=sf(1,1:N_m);
sfw=sf(1,N_m+1:end);

%%
% if v > sf, this function returns x = 1. Otherwise returns x = 0
x = 1;

for k = 1:N_m
    if sfm(1,k)~=0
        if vm(1,k)==0
            x = 0;
            break
        elseif M(vm(1,k),k) < M(sfm(1,k),k)
            x = 0;
            break
        end
    end
end

if x == 1
    for k = 1:N_w
        if sfw(1,k) ~= 0
            if vw(1,k) == 0
                x = 0;
                break
            elseif W(vw(1,k),k)<W(sfw(1,k),k)
                x = 0;
                break
            end
        end
    end
end

if x == 1
    if isequal(v, sf) == 1
        x = 0;
    end
end




