function [y] = f(t,f0)
%generating the function in exercise B
T = 1/f0;
y = zeros(1,length(t));

%implementing the definition of the function
for time_index = 1:length(t)
    if mod(t(time_index),T) < T/2
        y(time_index) = 2/T*t(time_index);
    else
        y(time_index) = 0;
    end
end