function H = RComplementHamming(x,y)
% 计算整数向量x与y的反补汉明距离,返回整数值d
% x:整数行向量
% y:整数行向量
% x与y必须相同长度
y = fliplr(3-y);
H = size(find(abs(x-y)>0),2);
end