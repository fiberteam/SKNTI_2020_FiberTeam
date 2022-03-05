function [outputArg1,outputArg2] = MyEye(x,s,samplerate)
x=x(:)';
k = floor(length(x)/s);

t = linspace(-1, 1, s);
t = t / samplerate;

for i=5:k
    a = (i-1)*s+1;
    b = a + s -1;
    plot(t, x(a:b), 'color', 'blue');
end
xlim([t(1) t(end)])
xlabel('t [s]')
ylabel('amplitude')


end

