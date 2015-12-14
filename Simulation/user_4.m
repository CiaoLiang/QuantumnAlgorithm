clear all
clc
%%
monte_carlo_number = 100;% Number of Monte Carlo
K = 4;% Number of user
x = [0 0 0 0;0 0 0 1;0 0 1 0;0 0 1 1;0 1 0 0;0 1 0 1;0 1 1 0;0 1 1 1
     1 0 0 0;1 0 0 1;1 0 1 0;1 0 1 1;1 1 0 0;1 1 0 1;1 1 1 0;1 1 1 1
    ];
L_info = 2^11;% infomation bits per frame
L_tail = 0;
L_total = L_info + L_tail;% total bits per frame
niter = 20;% Number of iterations
c_length = 4 ;% spreading length
c = repmat([1,-1],1,c_length/2);% Spreading sequence
Rc = 1/c_length;% The rate
n_effi = 1.0;%photoelectric conversion efficiency
nb_s = 39/n_effi;%average background photons for each bit
ns_b = 40;   %average signal photons for each bit
nb = nb_s*Rc;%average background photons for each chip
M1 = 2*ns_b*Rc;%photons count for chip 1
M0 = 0;% photon counts for chip 0
h = ones(1,K);%default Channel coefficient
%%
fprintf('\n\n=====================IDMA    MUD    Demo======================\n');
fprintf(' Whole Block size = %6d\n',L_total);
fprintf(' Spreading sequence length = %6d\n',c_length);
fprintf(' photoelectric conversion efficiency = %6d\n',n_effi);
fprintf(' iteration number =  %6d\n', niter);
fprintf(' monte carlo number = %6d\n', monte_carlo_number);
fprintf(' equivalent background radiation photons = %6d\n', nb_s);
fprintf(' user number =  %6d\n', K);
fprintf(' ns = ');
fprintf('\n');
for i = 1:length(ns_b)
    fprintf('%15.2f ',ns_b(i));
    if mod(i,4)==0
        fprintf('\n');
    end
end
if mod(length(ns_b),4)~=0
    fprintf('\n');
end
fprintf('===============================================================\n');
fprintf('BER = \n')
fprintf('\n===============================================================\n');
fprintf('>>>>>>>>>> Please be patient. Wait a while to get the result. >>>>>>>>>>\n');
%%
for i = 1:K
    [~, alpha1] = sort(rand(1,L_total*c_length));
    alpha(i,:)=alpha1;
end
clear  alpha1
%%
for nM1 = 1:length(M1)
    errs(:,:,nM1) = zeros(K,niter);
    nferr(:,:,nM1) = zeros(1,niter);
    nframe = 0;
    while nframe < monte_carlo_number
        nframe = nframe + 1
        d = round(rand(K, L_total));
        en_output=2*d-1;% channel encoder
        spread_en_output = (signal_spread( en_output,c )+1)/2;% spread
        clear en_output;
        for k = 1:K
            inter_spread_en_output(k,:) = spread_en_output(k,alpha(k,:));%interleave
            for i=1:L_total*c_length
                if inter_spread_en_output(k,i)==0
                    inter_spread_en_output(k,i)=M0;
                else
                    inter_spread_en_output(k,i)=M1(nM1);
                end
            end
        end;        
        clear spread_en_output;
        r = poissrnd(n_effi*sum(diag(h) * inter_spread_en_output,1)+nb);
%%
        L_a = zeros( K, L_total*c_length);% Initialize a priori LLR
        for iter = 1:niter
          %% MUD
            E_s=exp(L_a)./(1+exp(L_a));
            est_noise=n_effi * M1(nM1) * repmat(sum(diag(h)*E_s,1),K,1)+ nb-n_effi * M1(nM1)*diag(h)*E_s;
            lamda_1=est_noise+n_effi*M1(nM1)*repmat(h',1,L_total*c_length);
            lamda_0=est_noise;
            fai_1=repmat(r,K,1).*log(lamda_1)-lamda_1;
            fai_0=repmat(r,K,1).*log(lamda_0)-lamda_0;
            chi_1=zeros(K,L_total*c_length);
            chi_0=zeros(K,L_total*c_length);
            sum_1=zeros(K,L_total*c_length);
            sum_0=zeros(K,L_total*c_length);
            for t = 1 : K^2
                for k = 1 : K
                    chi_1(k,:)=exp(sum(repmat(x(t,:)',1,L_total*c_length).*lamda_1,1)...
                        -repmat(x(t,k),1,L_total*c_length).*L_a(k,:));
                    chi_0(k,:)=exp(sum(repmat(x(t,:)',1,L_total*c_length).*lamda_1,1)...
                        -repmat(x(t,k),1,L_total*c_length).*L_a(k,:));
                end;
                sum_1=sum_1+chi_1;
                sum_0=sum_0+chi_0;
            end;
            sum_1=sum_1/2;
            sum_0=sum_0/2;
            L_e = log(sum_1)-log(sum_0)+fai_1-fai_0;
          %%
            for k = 1:K
                L_a(k,alpha(k,:)) = L_e(k,:);
                L_a_1=round((sign(L_a-700)+1)/2);
                L_a=700*L_a_1+L_a.*(1-L_a_1);
            end
          %% DEC
            for k = 1:K
                [L_e(k,:),estimate_d_k(k,:)] = idma_app( L_a(k,:),c);
            end
          %%
            for k = 1:K
                L_a(k,:) = L_e(k,alpha(k,:));
            end
          %%
            est_d = estimate_d_k(:,1:L_info);% ignore tail bits
            %Counterr bit errors for the current iteration
            err(:,iter) = zeros(K,1);
            for k = 1:K
                err(k,iter) = err(k,iter) + length(find(est_d(k,:)~=d(k,:)));
            end
            %Count frame errors for the current iteration
            if sum(err(:,iter))>0
                nferr(1,iter,nM1) = nferr(1,iter,nM1)+1;
            end
        end;%for iter    
%%
        for k = 1:K
            errs(k,1:niter,nM1) = errs(k,1:niter,nM1) + err(k,1:niter);
        end
        if nframe==monte_carlo_number
            ber(1:K,1:niter,nM1) = errs(1:K,1:niter,nM1)/(nframe*L_info);% Bit error rate
            total_ber(1,1:niter,nM1) = sum(errs(1:K,1:niter,nM1),1)/(nframe*K*L_info);% Total Bit error rate
            fer(1,1:niter,nM1) = nferr(1,1:niter,nM1)/nframe; % Frame error rate
            fprintf('==================== ns = %d ============================\n', ns_b(nM1))
            fprintf('%d frames transmitted, %d frames in error.\n', nframe, nferr(1,niter,nM1));
            fprintf('%d bits were transmitted totally.\n',nframe*K*L_info);
            fprintf('Bit Error Rate (from iteration 1 to iteration %d):\n', niter);
            for i=1:niter
                fprintf('%8.4e    ', total_ber(1,i,nM1));
                if mod(i,5)==0
                    fprintf('\n');
                end
            end
            if mod(niter,5)~=0
                fprintf('\n');
            end
        end
    end%while
end%nM1
figure;
title('Operation is completed, please check the results.');




