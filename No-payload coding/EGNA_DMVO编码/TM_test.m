

A=newDNA;
B=char(A+'0')
[m,l]=size(B);

for i=1:m
    C(i,:)=strrep(B(i,:),'0','T') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'1','C') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'2','G') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'3','A') ;
end
[Tm,GC]=GCTmBioBox(C)
%[Tm1,GC1]=GCTmBioBox(newDNA);

std(Tm,0)
%std(Tm1,0)
clear all 
