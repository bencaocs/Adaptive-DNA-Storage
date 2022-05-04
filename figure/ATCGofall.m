clear;
% %ֻ����ʽ��txt�ļ�
% file_t = fopen('mytxt.txt','r');
% %��ʮ���ƶ�ȡ,�Ҷ�ȡ�������Զ��ų�һ�У��ŵ�˳��Ϊ���ȴӵ�һ����ߵ���һ���ұߣ�Ȼ���ŵڶ���
% A = fscanf(file_t,'%d');
% %�ر��ļ�
% fclose(file_t);
data=importdata('Encoded_oligo_60_0.45-0.55_16050_18000.txt');
data1=importdata('Encoded_oligo_16050_18000.txt');
%ATCG=['TAAGAGCCACTGGGAACAGAACGTGGATAATCCTTAGTGTTCGGGTAGATGAAGGACGCTCTACGCCGTTAAGACGGCATGATATACAAAGTCATACAGACGTTTCTACCACGGCTAGCCTAAGGGTCAAGAGGATAAGAAGACTTTGATCC'];

%Draw the search space   [Y,I]=sort(a,'descend')
getall_pic1 = juzhenfenbu(data);
getall_pic2 = juzhenfenbu(data1);
%allpic=[sort(abs(getall_pic1(1:5000,1)),'descend'),sort(abs(getall_pic2(1:5000,1)),'descend')];
allpic=[getall_pic1(800:1300,3) getall_pic2(800:1300,3)];
plot(allpic);
legend('Balance bases','old')
plot(getall_pic1(1:200,:));
plot(getall_pic1(1:100,2));
title('Encoded A/T/G/C-content')
xlabel('Article sequence number');
ylabel('A/T/G/C-concent');



axis tight
grid off
box on
legend('base-concent')

