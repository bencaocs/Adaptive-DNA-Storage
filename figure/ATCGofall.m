clear;
% %只读形式打开txt文件
% file_t = fopen('mytxt.txt','r');
% %以十进制读取,且读取的数据自动排成一列，排的顺序为：先从第一行左边到第一行右边，然后排第二行
% A = fscanf(file_t,'%d');
% %关闭文件
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

