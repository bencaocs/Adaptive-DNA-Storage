function EditDistMatrix= EditDistMatrix(IntDNAMatrix,R)
%�������������DNA���������м�ĺ�������,���ؾ������(�Գ�)
%IntDNAMatrix:��������,��������DNA���о���
%R:�ṹ��,�������
%DistHammingMatrix:��������(�Գ�),a_ij=HammingDist(DNA_i,DNA_j)
%HammingDistMatrix = pdist2(IntDNAMatrix,IntDNAMatrix,'hamming')*size(IntDNAMatrix,2);
%HammingDistMatrix =EditDist(IntDNAMatrix,IntDNAMatrix);
for i=1:size(IntDNAMatrix,1)
    for j=1:size(IntDNAMatrix,1)
    EditDistMatrix(i,j)=EditDist(IntDNAMatrix(i,:),IntDNAMatrix(j,:));
    end
end
end