function EditDistcon= Editconstraint(d,DNAset)
%�������������DNA���������м�ı༭����,����һ������
%IntDNAMatrix:��������,��������DNA���о���
%R:�ṹ��,�������
    len=size(d,2);
    a=len-2;
    b=len-1;
%IntDNAMatrix = cat(1,d,IntDNAMatrix);
% Hamming = pdist2(IntDNAMatrix,IntDNAMatrix,'hamming')*size(IntDNAMatrix,2);
% HammingSUM = sum(Hamming(1,:));
for i=1:size(DNAset,1)

    EditDistcon(1,i)=EditDist(DNAset(i,:),d);
    EditFrist(1,i)=EditDist(DNAset(i,[1,2,3]),d(len-2:len-1:len));
    EditLast(1,i)=EditDist(DNAset(i,[a b len]),d([1 2 3]));
    if( EditLast(1,i)==0|| EditFrist(1,i)==0)
     EditDistcon(1,i)=0;
    end
    
end

end