function H = RCHamming(x,IntDNAMatrix)
% ������������x��y�ĺ�������,��������ֵd
% x:����������
% y:����������
% x��y������ͬ����
IntDNAMatrix = fliplr(3-IntDNAMatrix);
 H = pdist2(x,IntDNAMatrix,'hamming')*size(IntDNAMatrix,2);
end