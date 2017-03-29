function [fselect, stftselect, L, P, M, I] = get_stft_features(F, S, tH)
select = F>tH(1) & F<tH(2);
fselect = F(select);
stftselect = S(select);
[M, I] = max(stftselect);
if length(fselect)>2
    [P1, L1] = findpeaks(stftselect,'SortStr','descend');
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