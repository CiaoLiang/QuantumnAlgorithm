diary idma_mud_demo_bpsk2.txt
clear all
%=====================K=8=============1500 frames==========================
L_info = 256;
L_tail = 0;
L_total = L_info + L_tail;
niter = 4;%修改作为不同迭代次数，得到仿真性能进行比较
monte_carlo_number = 100;
EbN0db = [4];%可设置为数组，比较不同信噪比下的仿真性能比较
K = 3;%用户数要大一点，逼真度更好，但须较大信噪比支持
c = repmat( [1,-1], 1, 5);
c_length = length(c);
Rc = 1/c_length;
flag = 0; %1 Orignal 0 QuantumAssistant

% h = ones(1,K);  % AWGN channel信道参数h是1？？
% h =  randn(1,K);  % Rayleigh Channel

for i = 1:K
    [temp, alpha1] = sort(rand(1,L_total*c_length));   %将随机向量元素按照升序排列成temp；  alphal是temp的索引，
    alpha(i,:)=alpha1;                                 %每次循环生成一个交织器，共K个
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
for i = 1:length(EbN0db)   %这里是输出 信噪比数组值
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
    sigma = 1/sqrt(2*Rc*en); 	   % standard deviation of AWGN noise       @@@@@@@@@@@是这样计算
    
  % Clear bit error counter and frame error counter
    errs(:,:,nEN) = zeros(K,niter);%在while循环之外就是为了每次循环不清除原有的统计值，5帧一统计，叠加统计
    nferr(:,:,nEN) = zeros(1,niter);
        
  % clear counter of transmitted frames  
    nframe = 0;    
    while nframe < monte_carlo_number
          nframe = nframe + 1;
          d = round(rand(K, L_info));    % info. bits，生成用户信息比特,rand产生均匀分布的随机数组（0和1范围内），round是四舍五入至整数
          en_output = d ; % channel encoder
          en_output = 2*en_output-1;     % symbol mapper,将0转变为-1，将1转变为+1
          spread_en_output = signal_spread( en_output,c );    % spread ，   调用扩频函数         
          
          % h = ones(1,K);%信道参数全为1
           h =  randn(1,K);%rayley信道
          % h = ones(1,K);
          for k = 1:K
              inter_spread_en_output(k,:) = spread_en_output(k,alpha(k,:)); % interleaver，将每个用户的扩频后码片分别进行交织处理
          end
          r = sum(diag(h) * inter_spread_en_output,1) + sigma * randn(1,L_total*c_length) ; % received bits,将交织后所有用户序列乘以信道系数h按列加起来形成一行，然后加上模拟的加性高斯白噪声，之和作为接收序列
          fn = zeros(2^K,L_total*c_length); %fn是最原始的CostFunction
          
          for i = 1:2^K
                index = IndexToBinary(i-1,K);
                fn(i,:) = exp((-(r-h*repmat(index',1,L_total*c_length)).^2)/(2*sigma^2));               
          end %计算不带f(x)的CostFunction
       
          % Initialize extrinsic information      
          L_e = zeros( K, L_total*c_length); %这里的L_e即先验信息  L ESE()
          
          for iter = 1:niter
              % ESE operations
              L_a = L_e;  % a priori info. L_a是先验信息,也是 L ESE()
  
              if flag == 0 %QuanntumAlogrithm
                  P_0 = 1./(1+exp(L_a));
                  P_1 = 1-P_0; %为0和1的概率,用于计算P(X),X是向量矩阵
               P_0 = P_0*10;
               P_1 = P_1*10; %计算提高精度，中间可能存在精度的误差
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
              end %计算CostFunction

              X_0 = zeros(2^K/2,L_total*c_length,K); %0解集
              X_1 = zeros(2^K/2,L_total*c_length,K); %1解集
              
                  m0 = ones(1,K); % 计数标记
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
                  end %生成子解集

                  L_po = zeros( K,L_total*c_length);%初始化L_po
                  for k = 1:K
                     L_po(k,:) = log((sum(X_1(:,:,k)))./(sum(X_0(:,:,k)))); %S0-DHA QMUD With Maximum approximation without DHA
                  end
                  %注意别比反了,X是指发送X向量的概率。这个地方要改。
                  %初步采用sum方法
                  L_e = L_po-L_a;
              end 
              
              if flag == 1 %Orignal
                  E_x = tanh(L_a/2);
                  Var_x = 1-E_x.^2;
                  E_r = sum(diag(h)*E_x,1);
                  Var_r = sum(diag(h.^2)*Var_x,1)+sigma^2;
                  L_e = 2*diag(h)*(repmat(r,K,1)-repmat(E_r,K,1)+diag(h)*E_x)./(repmat(Var_r,K,1)-diag(h.^2)*Var_x);%这里的L_e指的是外信息e ESE()!!!!      repmat使行向量按行扩展，从而每个用户都使用同样的一个接收序列
              end 
                  
              % a priori info. and de-interleave
              for k = 1:K
                  L_a(k,alpha(k,:)) = L_e(k,:);     % e ESE()经过解交织，成为DEC的输入先验信息 L DEC(),即这里的L_a指的是L DEC()
              end
             
              % DEC operations
              for k = 1:K
                  [L_SISO_2(k,:),estimate_d_k(k,:)] = idma_app( L_a(k,:),c); %DEC 标准的APP解码，调用解码函数，L_SISO_2指的是e DEC()，estimate_d_k指的是译码输出序列
              end
              
              % interleave
              for k = 1:K
                  L_e(k,:) = L_SISO_2(k,alpha(k,:));   %e DEC()经过交织，变为L ESE()，这里的L_e即  L ESE()，从而完成一次迭代
              end 
              
              % Number of bit errors in current iteration
              est_d = estimate_d_k(:,1:L_info); % ignore tail bits
              err(:,iter) = zeros(K,1);
              for k = 1:K
                  err(k,iter) = err(k,iter) + length(find(est_d(k,:)~=d(k,:)));%find返回索引向量。完成第iter次迭代中所有K用户各自的误比特数目
              end
              % Count frame errors for the current iteration
              if sum(err(:,iter))>0
                  nferr(1,iter,nEN) = nferr(1,iter,nEN)+1;%所有K用户的信息比特都没有错误传输时才视为 一帧无差错传输，否则统计错误帧数
              end   
          end   %iter
          
          % Total number of bit errors for all iterations
          for k = 1:K
              errs(k,1:niter,nEN) = errs(k,1:niter,nEN) + err(k,1:niter);%实际上将每一帧处理得到的errs数目累加起来，便于后面与帧数nframe相除得到多帧的统计平均
          end
          
          if rem(nframe,1000)==0 | nframe==monte_carlo_number       %只有当传输的帧的数目是5的倍数或蒙特卡洛时，才做以下操作，即每隔5帧一输出统计信息
              % Bit error rate
              ber(1:K,1:niter,nEN) = errs(1:K,1:niter,nEN)/(nframe*L_info);%得到每次迭代每个用户的比特误码率
              % Total Bit error rate
              total_ber(1,1:niter,nEN) = sum(errs(1:K,1:niter,nEN),1)/(nframe*K*L_info);%每次迭代所有用户总的平均比特误码率
              % Frame error rate
              fer(1,1:niter,nEN) = nferr(1,1:niter,nEN)/nframe;%每次迭代帧错误率
          
              % Display intermediate results in process
              fprintf('******************* Eb/N0 = %2.1f dB *********** K = %2d ********************\n', EbN0db(nEN),K);
              fprintf('%d frames transmitted, %d frames in error.\n', nframe, nferr(1,niter,nEN));%最后一次迭代得到的帧错误数！！！！
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
                  fprintf('%8.4e    ', total_ber(1,i,nEN));%每次迭代所有用户总的平均比特误码率
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