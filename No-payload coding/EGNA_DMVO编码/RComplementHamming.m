function H = RComplementHamming(x,y)
% ������������x��y�ķ�����������,��������ֵd
% x:����������
% y:����������
% x��y������ͬ����
y = fliplr(3-y);
H = size(find(abs(x-y)>0),2);
end