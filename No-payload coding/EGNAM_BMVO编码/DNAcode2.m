function np=DNAcode2(p)
[m,l]=size(p);
for i=1:m
    for j=1:l %j��ֵ�ǣ���ʼֵ��1,��2Ϊ��������py����
        switch p(i,j) 
            case 0
                np(i,j)= 'T';
            case 1
                np(i,j)= 'C';
            case 2
                np(i,j)= 'G';
            case 3
                np(i,j)= 'A';
        end
    end
end

