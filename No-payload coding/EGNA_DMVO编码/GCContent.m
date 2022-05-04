function  flag = GCContent(candidate)
[px,py] = size(candidate);
count = 0;
for i = 1:py
    if candidate(1,i) == 1 || candidate(1,i) == 2
           count = count + 1;
    end
end
if count == floor(py/2)
    flag = 1;
else
    flag =0;
end
