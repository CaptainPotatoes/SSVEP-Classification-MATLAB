function [ count ] = countOccurrences( data, CLASS )
data = data(:); %Vectorize;
count = 0;
for i=1:length(data)
    if data(i) == CLASS
        count = count+1;
    end
end

end

