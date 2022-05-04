% function r = testResult(result,d)
% A=RCDisHammingMatrix(result);
% B=DistHammingMatrix(result);
% B(logical(eye(size(B))))=d;
% px=size(result,1);
% temp=[];
% for i=1:px
%  temp(i)=GC(result(i,:));
% end
%  if isempty(find(A<d, 1)) && isempty(find(B<d, 1)) &&  sum(temp)==px
%        r=1;
%  else
%      r = 0;
%  end
Param1=[' ATGCATCGATCAG'];%第一个参数，一定要有' '
Param2=[' TACGTAGCTAGTC']; 
MFE_fun(Param1,Param2)



