diary idma_mud_demo_bpsk2.txt
clear all
%=====================K=8=============1500 frames==========================
L_info = 256;
L_tail = 0;
L_total = L_info + L_tail;
niter = 4;%�޸���Ϊ��ͬ�����������õ��������ܽ��бȽ�
monte_carlo_number = 100;
EbN0db = [4];%������Ϊ���飬�Ƚϲ�ͬ������µķ������ܱȽ�
K = 3;%�û���Ҫ��һ�㣬����ȸ��ã�����ϴ������֧��
c = repmat( [1,-1], 1, 5);
c_length = length(c);
Rc = 1/c_length;
flag = 0; %1 Orignal 0 QuantumAssistant

% h = ones(1,K);  % AWGN channel�ŵ�����h��1����
% h =  randn(1,K);  % Rayleigh Channel

for i = 1:K
    [temp, alpha1] = sort(rand(1,L_total*c_length));   %���������Ԫ�ذ����������г�temp��  alphal��temp��������
    alpha(i,:)=alpha1;                                 %ÿ��ѭ������һ����֯������K��
end
clear temp alpha1
%==========================================================================

fprintf('\n\n=====================IDMA    MUD    Demo======================\n'); 
fprintf(' Whole Block size = %6d\n',L_total);
fprintf(' Spreading sequence length = %6d\n',c_length);
fprintf(' iteration number =  %6d\n', niter);
fprintf(' monte carlo number = %6d\n', monte_carlo_number);
fprintf(' Eb / N0 (dB) = ');
fprintf('\n');
for i = 1:length(EbN0db)   %��������� ���������ֵ
    fprintf('%10.2f ',EbN0db(i));
    if mod(i,4)==0
        fprintf('\n');
    end
end
if mod(length(EbN0db),4)~=0
    fprintf('\n');
end
fprintf(' user number =  %6d', K);
%fprintf(' Channel coefficient: \n');
%for i = 1:K
%    fprintf( '%6d', h(i));
%end
%fprintf('\n')
fprintf('\n===============================================================\n\n');

fprintf('>>>>>>>>>> Please be patient. Wait a while to get the result. >>>>>>>>>>\n');

for nEN = 1:length(EbN0db)
    en = 10^(EbN0db(nEN)/10);      % convert Eb/N0 from unit db to normal numbers
    sigma = 1/sqrt(2*Rc*en); 	   % standard deviation of AWGN noise       @@@@@@@@@@@����������
    
  % Clear bit error counter and frame error counter
    errs(:,:,nEN) = zeros(K,niter);%��whileѭ��֮�����Ϊ��ÿ��ѭ�������ԭ�е�ͳ��ֵ��5֡һͳ�ƣ�����ͳ��
    nferr(:,:,nEN) = zeros(1,niter);
        
  % clear counter of transmitted frames  
    nframe = 0;    
    while nframe < monte_carlo_number
          nframe = nframe + 1;
          d = round(rand(K, L_info));    % info. bits�������û���Ϣ����,rand�������ȷֲ���������飨0��1��Χ�ڣ���round����������������
          en_output = d ; % channel encoder
          en_output = 2*en_output-1;     % symbol mapper,��0ת��Ϊ-1����1ת��Ϊ+1
          spread_en_output = signal_spread( en_output,c );    % spread ��   ������Ƶ����         
          
          % h = ones(1,K);%�ŵ�����ȫΪ1
           h =  randn(1,K);%rayley�ŵ�
          % h = ones(1,K);
          for k = 1:K
              inter_spread_en_output(k,:) = spread_en_output(k,alpha(k,:)); % interleaver����ÿ���û�����Ƶ����Ƭ�ֱ���н�֯����
          end
          r = sum(diag(h) * inter_spread_en_output,1) + sigma * randn(1,L_total*c_length) ; % received bits,����֯�������û����г����ŵ�ϵ��h���м������γ�һ�У�Ȼ�����ģ��ļ��Ը�˹��������֮����Ϊ��������
          fn = zeros(2^K,L_total*c_length); %fn����ԭʼ��CostFunction
          
          for i = 1:2^K
                index = IndexToBinary(i-1,K);
                fn(i,:) = exp((-(r-h*repmat(index',1,L_total*c_length)).^2)/(2*sigma^2));               
          end %���㲻��f(x)��CostFunction
       
          % Initialize extrinsic information      
          L_e = zeros( K, L_total*c_length); %�����L_e��������Ϣ  L ESE()
          
          for iter = 1:niter
              % ESE operations
              L_a = L_e;  % a priori info. L_a��������Ϣ,Ҳ�� L ESE()
  
              if flag == 0 %QuanntumAlogrithm
                  P_0 = 1./(1+exp(L_a));
                  P_1 = 1-P_0; %Ϊ0��1�ĸ���,���ڼ���P(X),X����������
               P_0 = P_0*10;
               P_1 = P_1*10; %������߾��ȣ��м���ܴ��ھ��ȵ����
              cf = zeros(2^K,L_total*c_length);
              Pro_X = ones(2^K,L_total*c_length);

              for i = 1:2^K
                index = IndexToBinary(i-1,K);
                for j = 1:K
                    for k = 1:L_total*c_length
                        if index(j) == -1
                            Pro_X(i,k) = Pro_X(i,k).*P_0(j,k);
                        end
                        if index(j) == 1
                            Pro_X(i,k) = Pro_X(i,k).*P_1(j,k);
                        end
                    end
                end
                cf(i,:) = fn(i,:).*Pro_X(i,:);               
              end %����CostFunction

              X_0 = zeros(2^K/2,L_total*c_length,K); %0�⼯
              X_1 = zeros(2^K/2,L_total*c_length,K); %1�⼯
              
                  m0 = ones(1,K); % �������
                  m1 = ones(1,K);
                  for i = 1:2^K         
                      index = IndexToBinary(i-1,K);
                      for k = 1:K
                          if index(k) == -1
                               X_0(m0(k),:,k) = cf(i,:);
                               m0(k)=m0(k)+1;
                          end
                          if index(k) == 1
                               X_1(m1(k),:,k) = cf(i,:);
                               m1(k)=m1(k)+1;
                          end
                      end
                  end %�����ӽ⼯

                  L_po = zeros( K,L_total*c_length);%��ʼ��L_po
                  for k = 1:K
                     L_po(k,:) = log((sum(X_1(:,:,k)))./(sum(X_0(:,:,k)))); %S0-DHA QMUD With Maximum approximation without DHA
                  end
                  %ע���ȷ���,X��ָ����X�����ĸ��ʡ�����ط�Ҫ�ġ�
                  %��������sum����
                  L_e = L_po-L_a;
              end 
              
              if flag == 1 %Orignal
                  E_x = tanh(L_a/2);
                  Var_x = 1-E_x.^2;
                  E_r = sum(diag(h)*E_x,1);
                  Var_r = sum(diag(h.^2)*Var_x,1)+sigma^2;
                  L_e = 2*diag(h)*(repmat(r,K,1)-repmat(E_r,K,1)+diag(h)*E_x)./(repmat(Var_r,K,1)-diag(h.^2)*Var_x);%�����L_eָ��������Ϣe ESE()!!!!      repmatʹ������������չ���Ӷ�ÿ���û���ʹ��ͬ����һ����������
              end 
                  
              % a priori info. and de-interleave
              for k = 1:K
                  L_a(k,alpha(k,:)) = L_e(k,:);     % e ESE()�����⽻֯����ΪDEC������������Ϣ L DEC(),�������L_aָ����L DEC()
              end
             
              % DEC operations
              for k = 1:K
                  [L_SISO_2(k,:),estimate_d_k(k,:)] = idma_app( L_a(k,:),c); %DEC ��׼��APP���룬���ý��뺯����L_SISO_2ָ����e DEC()��estimate_d_kָ���������������
              end
              
              % interleave
              for k = 1:K
                  L_e(k,:) = L_SISO_2(k,alpha(k,:));   %e DEC()������֯����ΪL ESE()�������L_e��  L ESE()���Ӷ����һ�ε���
              end 
              
              % Number of bit errors in current iteration
              est_d = estimate_d_k(:,1:L_info); % ignore tail bits
              err(:,iter) = zeros(K,1);
              for k = 1:K
                  err(k,iter) = err(k,iter) + length(find(est_d(k,:)~=d(k,:)));%find����������������ɵ�iter�ε���������K�û����Ե��������Ŀ
              end
              % Count frame errors for the current iteration
              if sum(err(:,iter))>0
                  nferr(1,iter,nEN) = nferr(1,iter,nEN)+1;%����K�û�����Ϣ���ض�û�д�����ʱ����Ϊ һ֡�޲���䣬����ͳ�ƴ���֡��
              end   
          end   %iter
          
          % Total number of bit errors for all iterations
          for k = 1:K
              errs(k,1:niter,nEN) = errs(k,1:niter,nEN) + err(k,1:niter);%ʵ���Ͻ�ÿһ֡����õ���errs��Ŀ�ۼ����������ں�����֡��nframe����õ���֡��ͳ��ƽ��
          end
          
          if rem(nframe,1000)==0 | nframe==monte_carlo_number       %ֻ�е������֡����Ŀ��5�ı��������ؿ���ʱ���������²�������ÿ��5֡һ���ͳ����Ϣ
              % Bit error rate
              ber(1:K,1:niter,nEN) = errs(1:K,1:niter,nEN)/(nframe*L_info);%�õ�ÿ�ε���ÿ���û��ı���������
              % Total Bit error rate
              total_ber(1,1:niter,nEN) = sum(errs(1:K,1:niter,nEN),1)/(nframe*K*L_info);%ÿ�ε��������û��ܵ�ƽ������������
              % Frame error rate
              fer(1,1:niter,nEN) = nferr(1,1:niter,nEN)/nframe;%ÿ�ε���֡������
          
              % Display intermediate results in process
              fprintf('******************* Eb/N0 = %2.1f dB *********** K = %2d ********************\n', EbN0db(nEN),K);
              fprintf('%d frames transmitted, %d frames in error.\n', nframe, nferr(1,niter,nEN));%���һ�ε����õ���֡��������������
              %for k = 1:K
              %    fprintf('User %d Bit Error Rate (from iteration 1 to iteration %d):\n', k,niter);
              %    for i=1:niter
              %        fprintf('%8.4e    ', ber(k,i,nEN));
              %        if mod(i,5)==0
              %            fprintf('\n');
              %        end
              %    end
              %    if mod(niter,5)~=0
              %        fprintf('\n');
              %    end
              %end
              
              fprintf('Bit Error Rate (from iteration 1 to iteration %d):\n', niter);
              for i=1:niter
                  fprintf('%8.4e    ', total_ber(1,i,nEN));%ÿ�ε��������û��ܵ�ƽ������������
                  if mod(i,5)==0
                      fprintf('\n');
                  end
              end
              if mod(niter,5)~=0
                  fprintf('\n');
              end
              
              %fprintf('Frame Error Rate (from iteration 1 to iteration %d):\n', niter);
              %for i=1:niter
              %    fprintf('%8.4e    ', fer(1,i,nEN));
              %    if mod(i,5)==0
              %        fprintf('\n');
              %    end
              %end
              %if mod(niter,5)~=0
              %    fprintf('\n');
              %end
              
              fprintf('**************************************************************************\n\n');
              save idma_mud_demo_bpsk2 EbN0db ber fer total_ber
          end  %if
      end		%while
end 		%nEN
diary off