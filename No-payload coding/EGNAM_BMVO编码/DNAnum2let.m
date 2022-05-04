function C = DNAnum2let(seq1)
B=char(seq1+'0');
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
i=1;

for i=1:m
    C(i,:)=strrep(C(i,:),'T','3') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'C','2') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'G','1') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'A','0') ;
end
i=1;
for i=1:m
    C(i,:)=strrep(C(i,:),'0','T') ;
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

    C=fliplr(C);

end