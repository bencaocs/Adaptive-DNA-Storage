function p1 =repair(popsize,chromlength,p1)
for i=1:popsize
    tempp1 = find(p1(i,1:chromlength) == 2 | p1(i,1:chromlength) == 1);
    tempp2 = find(p1(i,1:chromlength) == 0 | p1(i,1:chromlength) == 3);
    templen1=length(tempp1);
    templen2=length(tempp2);
    len = chromlength;
    repairlen = templen1 - floor(1/2*len);
    if repairlen ==0
        continue;
    end
    if repairlen > 0
        temparray = randperm(templen1);
        for j = 1:repairlen
            if rand(1) > 0.5
                p1(i,tempp1(temparray(j))) =0;
            else
                p1(i,tempp1(temparray(j))) =3;
            end
        end
    else 
        temparray = randperm(templen2);
        for j = 1:abs(repairlen)
            if rand(1) > 0.5
                p1(i,tempp2(temparray(j))) =1;
            else
                p1(i,tempp2(temparray(j))) =2;
            end
        end
    end
end

