%clear all    
     for i=1:size(newDNA,1)
         MFEvalue(i,:)=MFE_fun(DNAComplement(newDNA(i,:)),DNAnum2let(newDNA(i,:)));
     end
    
     MFEAve=sum(MFEvalue)/size(newDNA,1)
%           px = size(newDNA,1);
%      index2 = 1;
%      temp=[];
%      for q = 1:px
%          %if  MFE_fun(DNAComplement(newDNA(q,:)),DNAnum2let(newDNA(q,:)))>=MFEAve%||  GCContent(DNASet(q,:))~=1
%          if  MFEvalue(q,:)>=MFEAve%||  GCContent(DNASet(q,:))~=1
%              temp(index2) =q;
%              index2 = index2 + 1;
%          end
%      end
%      if length(temp)~=0
%          newDNA(temp,:)=[];
%      end
   
     fprintf('newDNA的大小为%d\n',size(newDNA,1));
     [Tm,GC]=GCTmBioBox(DNAComplement(newDNA));
     std(Tm,0) 
     clear all 