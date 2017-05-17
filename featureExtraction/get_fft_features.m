function [fselect, fftselect, L, P] = get_fft_features( f, FFT, threshold )
select = f>threshold(1) & f<threshold(2);
fselect = f(select);
fftselect = FFT(select);
[M, I] = max(fftselect);
if length(fselect)>2
    [P1, L1] = findpeaks(fftselect,'SortStr','descend');
else
    P1 = [];
    L1 = [];
end

if ~isempty(P1)
    L = fselect(L1(1));
    P = P1(1);
else
    L = 0;
    P = 0;
end

end