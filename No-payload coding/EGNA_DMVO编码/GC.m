function flag = GC(DNA)
[px,py] = size(DNA);
temp=[];
index = 1;
count=0;
for i=1:py
  if DNA(1,i) == 1 || DNA(1,i) == 2
           count = count + 1;
    end
end
if   count==floor(py/2)
    flag=1;
else
    flag=0;
end
 
end