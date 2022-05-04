function InitDNASet = NG_InitDNASet(Length,~,count)
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
    if   GCContent(DNASet(q,:))~=1
        temp(index) =q;
        index = index + 1;
    end
end
if length(temp)~=0
    DNASet(temp,:)=[];
end
%  Values= Continuity(DNASet(a,:));
%      a =1;
%      M = size(DNASet,1);
%   while a <= M
%     if Values > 1 
%            DNASet=[];
%            a = a-1;
%     end
%   end
% DNASet1 = Continuity(DNASet);%得到不能相同的个数       
% delete_index = [];
% for DNASet1_index = 1:size(DNASet1,1)
%     a = tabulate(DNASet(DNASet1_index,:)) ;
%     if DNASet1(DNASet1_index,:)~=0
%         if sum(a(:,2)>DNASet1(DNASet1_index,:))
%             delete_index(:,1) = DNASet1_index;
%         else
%         end
%     end
% end
%DNASet(delete_index,:) = [];
InitDNASet = DNASet;