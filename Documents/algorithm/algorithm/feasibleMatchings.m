function[cores] = feasibleMatchings(N_m, N_w, N_t)

%%N_m >= N_w

contracts = zeros(N_m,1+2*N_t);

for i=1:N_w
    for j=2:1+2*N_t
        contracts(i,j) = i + N_w*(j-2);
    end
end

combination=contracts(1,:);
for i=2:N_w
    combination=combvec(combination,contracts(i,:));
end

outsideOptions = zeros(N_m-N_w,size(combination,2));

combination = [combination ; outsideOptions];

combination = combination';


manside = zeros(0,N_m);
for i =1:size(combination,1)
    permutations = perms(combination(i,:));
    manside = [manside;permutations];
end

manside=unique(manside,'stable','rows');

womanside = zeros(size(manside,1),N_w);

for i = 1:size(manside,1)
    for j=1:N_m
        if manside(i,j)~=0
            x = manside(i,j);
            y = mod(x,N_w);
            if y == 0
                y = N_w;
            end
            z = ceil(x/N_w);
            % yさんはz番目のコントラクトでjとマッチしている。
            womanside(i,y)=(z-1)*N_m+j;
        end  
    end
    
end

matchings = [manside womanside];

cores = zeros(0,N_m+N_w);
for i=1:size(matchings,1)
    if stableOrNot_with_contract(matchings(i,:)) == 1
        cores = [cores; matchings(i,:)];
    end
end