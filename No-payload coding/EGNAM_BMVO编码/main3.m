clear all
Length =15;
d=3;
Universes_no=3000;

H_target=d;
count =1000;
Max_iteration=2000;

lb=0;
ub=3;
dim=Length;
%DNASet =  InitDNASet(Length,d,count);
newDNA = EGNA_BMVO(Universes_no,Max_iteration,lb,ub,dim,Length,H_target,count);%MVO(Universes_no,Max_iteration,lb,ub,dim);
% iter = 1;
% iter_max =1500;
% a=6e-5;
% b=1.5;
 
     for i=1:size(newDNA,1)
         MFEvalue(i,:)=MFE_fun(DNAComplement(newDNA(i,:)),DNAnum2let(newDNA(i,:)));
     end
    
     MFEAve=sum(MFEvalue)/size(newDNA,1);
    
     px = size(newDNA,1);
     index2 = 1;
     temp=[];
     for q = 1:px
         %if  MFE_fun(DNAComplement(newDNA(q,:)),DNAnum2let(newDNA(q,:)))>=MFEAve%||  GCContent(DNASet(q,:))~=1
         if  MFEvalue(q,:)>=MFEAve%||  GCContent(DNASet(q,:))~=1
             temp(index2) =q;
             index2 = index2 + 1;
         end
     end
     if length(temp)~=0
         newDNA(temp,:)=[];
     end
     fprintf('newDNA的大小为%d\n',size(newDNA,1));
         
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
