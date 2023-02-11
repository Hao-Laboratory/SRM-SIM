%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPE computes the deconvolved widefield image and iteratively    %
% estimates the structured pattern                                %
%                                                                 %
% Inputs:                                                         %
%         I_image_up  : upsampled measurements                    %
%         T_incoherent: incoherent transfer function              %
%         max_itr     : maximal iteration number                  %
%         reg_beta    : regularizer for widefield deconvolution   %
%         reg_delta   : regularizer for pattern estimation        %
%                                                                 %
% Outputs:                                                        %
%         Ip_est      : estimated patterns                        %
%         I_mean      : mean of the measurements (widefield)      %
%         I_mean_dec  : deconvolved widefield image               %
%                                                                 %
%           Copyright (C) 2016 Li-Hao Yeh                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ip_est,I_mean,I_mean_dec] = IPE(I_image_up,T_incoherent,max_itr,reg_beta,reg_delta)

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) real(fftshift(ifft2(ifftshift(x))));

[N,M,Nimg] = size(I_image_up);


I_mean = mean(I_image_up,3);
I_mean_dec = iF(F(I_mean).*T_incoherent./(T_incoherent.^2+reg_beta));
I_mean_dec(I_mean_dec<0) = 0;

I_mean_dec = I_mean_dec/max(I_mean_dec(:));

Ip_est = zeros(N,M,Nimg);

tic;
fprintf('| Pattern #  |  error ratio  | Elapsed time (sec) |\n');
for i = 1:Nimg
    forward = iF(T_incoherent.*F(Ip_est(:,:,i).*I_mean_dec));
    err_initial = sum(sum(abs(I_image_up(:,:,i) - forward).^2));
    for j = 1:max_itr
        forward = iF(T_incoherent.*F(Ip_est(:,:,i).*I_mean_dec));
        Ip_est(:,:,i) = Ip_est(:,:,i) + 2*I_mean_dec.*iF(T_incoherent.*F(I_image_up(:,:,i) - forward))/max(I_mean_dec(:))^2 ; % Gaussian gradient
        Ip_est(:,:,i) = iF(F(Ip_est(:,:,i)).*T_incoherent.^2./(T_incoherent.^2+reg_delta));
        temp = Ip_est(:,:,i);
        temp(temp<0) = 0;
        if j == 1            
            t = 1;
            Ip_est(:,:,i) = temp;
            tempp = temp;
        else
            tp = t;
            t = (1+sqrt(1+4*tp^2))/2;
            Ip_est(:,:,i) = temp + (tp-1)*(temp - tempp)/t;
            tempp = temp;
        end        
        
    end
    forward = iF(T_incoherent.*F(Ip_est(:,:,i).*I_mean_dec));
    err_final = sum(sum(abs(I_image_up(:,:,i) - forward).^2));
    fprintf('|    %2d      |     %.2f  %%   |        %.2f        |\n', i, err_final/err_initial*100,toc);
end
