function [y] = name_signal(t,f0)
%A function that returns a signal generated based on the name TUDOR
%T = 20, U = 21, D = 4, O = 15, R = 18

T = 1/f0;
y = zeros(1,length(t));

for time_index = 1:length(t)
    if mod(t(time_index),T) < T/5
        y(time_index) = 20;
    
    elseif mod(t(time_index),T) < 2*T/5 
        y(time_index) = 21;
    
    elseif mod(t(time_index),T) < 3*T/5 
        y(time_index) = 4;
    
    elseif mod(t(time_index),T) < 4*T/5 
        y(time_index) = 15;
    else
        y(time_index) = 18;
    end
end
