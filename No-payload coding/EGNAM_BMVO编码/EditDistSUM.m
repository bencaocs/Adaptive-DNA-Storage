function EditDistSUM= EditDistSUM(IntDNAMatrix,d)
%计算整数编码的DNA矩阵中序列间的汉明距离,返回距离矩阵(对称)
%IntDNAMatrix:整数矩阵,整数编码DNA序列矩阵
%R:结构体,冗余参数
%DistHammingMatrix:整数矩阵(对称),a_ij=HammingDist(DNA_i,DNA_j)
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