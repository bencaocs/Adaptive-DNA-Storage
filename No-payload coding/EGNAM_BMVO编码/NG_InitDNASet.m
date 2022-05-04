function [InitDNASet,MFEAve] = NG_InitDNASet(Length,~,count)
X_min = 0;
X_max = 3;
DNASet = X_min + round((X_max - X_min).*rand(count,Length));
%%保证相邻数字不一样
 for p = 1:size(DNASet,1)
   for q = 1:size(DNASet,2)-1
         if DNASet(p,q) == DNASet(p,q+1)
             del_idx(p,:) = p;
         end
    end
 end
 del_idx(del_idx(:,end)==0) = [];
 DNASet(del_idx,:) = [];
% DNASet =repair(count,len,DNASet);
% || GCContent(DNASet(i,:))~=s1
DNASet = unique(DNASet,'rows','stable');

px = size(DNASet,1);
index = 1;
temp=[];
for q = 1:px
    if   GCContent(DNASet(q,:))~=1 %|| MFEvalue(q,:)>=MFEAve
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
    DNASet(temp,:)=[];
end
%%保证MFE小于平均MFE
% for i=1:size(DNASet,1)
%     MFEvalue(i,:)=MFE_fun(fliplr(DNAComplement(DNASet(i,:))),DNAComplement(DNASet(i,:)));
% end
% MFEAve=sum(MFEvalue)/size(DNASet,1);
% px = size(DNASet,1);
% index = 1;
% temp=[];
% for q = 1:px
%     if  MFEvalue(q,:)>=MFEAve%||  GCContent(DNASet(q,:))~=1 
%         temp(index) =q;
%         index = index + 1;
%     end
% end
% if length(temp)~=0
%     DNASet(temp,:)=[];
% end
InitDNASet = DNASet;