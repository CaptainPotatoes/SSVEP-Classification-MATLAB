function [ Y ] = ssvepcfilt( X )
% [8 20] bandpass butterworth N=3
X=X(:);
b = [0.00259188624245653,0,-0.00777565872736959,0,0.00777565872736959,0,-0.00259188624245653];
a = [1,-5.12636114248572,11.2013241099204,-13.3475825890437,9.14852464832994,-3.42093110978683,0.545786834463780];
Y = filtfilt(b,a,X);
end

