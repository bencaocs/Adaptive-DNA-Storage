function EditDistSUM= EditDistSUM(IntDNAMatrix,d)
%�������������DNA���������м�ĺ�������,���ؾ������(�Գ�)
%IntDNAMatrix:��������,��������DNA���о���
%R:�ṹ��,�������
%DistHammingMatrix:��������(�Գ�),a_ij=HammingDist(DNA_i,DNA_j)
IntDNAMatrix = cat(1,d,IntDNAMatrix);
% Hamming = pdist2(IntDNAMatrix,IntDNAMatrix,'hamming')*size(IntDNAMatrix,2);
% HammingSUM = sum(Hamming(1,:));
for i=1:size(IntDNAMatrix,1)
    for j=1:size(IntDNAMatrix,1)
    EditDistMatrix(i,j)=EditDist(IntDNAMatrix(i,:),IntDNAMatrix(j,:));
    end
end
EditDistSUM = sum(EditDistMatrix(1,:));
end