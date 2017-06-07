function [ Y ] = convoluteSignals( X1, X2 )
% filters and convolutes 2 data channels. 
    fch = ssvepcfilt2(X1); %[5 40]
    fch2 = ssvepcfilt2(X2);
    conv2ch = conv(fch,fch2,'full');
    % if length is odd: 
    if  mod(length(conv2ch),2)==1
        Y = conv2ch(1:end-1);
    else
        Y = conv2ch;
    end
end

