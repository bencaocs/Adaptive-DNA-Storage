clear all
Length =8;
d=5;
Universes_no=1000;

H_target=d;
count =1000;
Max_iteration=1500;

lb=0;
ub=3;
dim=Length;
%DNASet =  InitDNASet(Length,d,count);
newDNA = EGNA_DMVO(Universes_no,Max_iteration,lb,ub,dim,Length,H_target,count);%MVO(Universes_no,Max_iteration,lb,ub,dim);
% iter = 1;
% iter_max =1500;
% a=6e-5;
% b=1.5;

         
DNASet=[];
while iter <= iter_max
    fprintf('进化第%d代,newDNA的大小为%d\n',iter,size(newDNA,1));
    pos = 1;
    DNASet =  InitDNASet(Length,d,10000+iter);
    len = size(DNASet,1);
    tempDNA=[];
    for i=1:len
        DisH=Hamming(DNASet(i,:),newDNA);
        len1=find(DisH<d);
       if isempty(len1) 
           newDNA=[DNASet(i,:);newDNA];
            elseif length(len1)==1  
                newDNA(len1,:)=[];
                newDNA=[DNASet(i,:);newDNA];
            end
    end
    iter= iter+1;
end
set = newDNA;
BestSet = DNAcode2(newDNA);
