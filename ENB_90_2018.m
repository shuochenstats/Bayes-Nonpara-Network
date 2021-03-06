%% Parameter initialization
G=90;
R=1:G;
para=[];
M=5000 % the length of the chain
K=ones(1,90) ; % initial networks
%K=[ones(1,30),ones(1,30)*2,ones(1,30)*3];
%K=1:90 

 
 

 
%% start MCMC...

W=X*X';

subj_n=30;
p=0.5; % p value is estimated by using local fdr
dgW=W(1,1)*(1/(1-p));% element subject to the size of W sub matrix 


for m=1:M, %mcmc loop, for each mcmc we record 2 things, Adjacency of GxG and lik all data...
    for g=1:G,
        %remove g from K, then assign to any existing cluster or new cluster and
        %calculate probability...
        Lik=[];K1=K;K1(g)=[]; %K1=the cluster remove g;
        K_label=unique(K1); % the labels of cluster;
        for c=1:length(unique(K1)), % assign region to existing clusters, and calculate likelihood.
            
            K(g)=K_label(c); %label the K...; we reset K here if g alone, the #K -1
            %assigned region change the cov of edges then calculate lik
            K_sets=[]; % how many edges are in cluster after assigning...
            nK=length(unique(K));
            Det=[];TR=[]; % reuse for each cluster assignment, a vector for a assignment
            for k=1:nK,  %for each cluster after assignment...
                set_k=find(K==K_label(k));% for each subset, we check how many edges are included
                %tic
                ind_pool=[]; % search index of the edges in the same cluster
                for i=1:length(I),

                    w=prod(single(ismember([I(i),J(i)],set_k)));
                    if w
                        ind_pool=[ind_pool,i]; % for example, 50 nodes will include 1225 edges in the communities
                    end

                end
                %toc
                bn1=length(ind_pool); %bn1=3

                % for likelihood and we need calculate the determint and
                % trace of tr(inv(Sigma)W) for the ind_pool...

                % log determinat:
                %math for determinat: \det(A)=\left(1-\rho \right)^{n-1}\left(1 + \rho n-\rho\right)
                %https://math.stackexchange.com/questions/226539/determinant-of-a-n-symmetric-square-matrix-with-diagonal-1
                % negative from inverse
                Det(k)=-subj_n/2*(log(1-p+bn1*p)+(bn1-1)*log(1-p)); % tested det 
                TR(k)=-.5*(dgW*bn1-(sum(sum(W(ind_pool,ind_pool))))*(p/(1-p)/(1-p+bn1*p)));
                % (sum(sum(W(ind_pool,ind_pool))))
                %TR(k)=-trace((1/(1-p))*(eye(bn1)-p/(1-p+bn1*p))*W(ind_pool,ind_pool))*0.5;
       
                %toc (1/(1-p)).*(eye(bn1)-p/(1-p+bn1*p))
                %inv((1-p).*eye(bn1)+p)
                K_sets=[K_sets, ind_pool];

            end % each luster
            % if g is assigned to a cluster, then the edges between clusters...
            ind_btw=find(ismember(1:length(I),K_sets)==0);
            lik_btw=-0.5*length(ind_btw)*W(1,1);
            Lik(c)=sum(Det)+sum(TR)+lik_btw;
            % measure the size of each edge cluster, computation. blk size
            %bn=cellfun(@length, K_sets);

        end % end assigning to each cluster...
        
        % Assigning to a new cluster...
        K(g)=max(K_label)+1;
        K_label=unique(K); %number of K_label increase 1.
        
        K_sets=[];
        nK=length(unique(K));
        Det=[];TR=[];
        for k=1:nK,  %for each cluster after assignment...
            set_k=find(K==K_label(k));% for each subset, we check how many edges are included
            %tic
            ind_pool=[]; % search index of the edges in the same cluster
            for i=1:length(I),

                w=prod(single(ismember([I(i),J(i)],set_k)));
                if w
                    ind_pool=[ind_pool,i]; % for example, 50 nodes will include 1225 edges in the communities
                end

            end
            bn1=length(ind_pool); %bn1=3

            % for likelihood and we need calculate the determint and
            % trace of tr(inv(Sigma)W) for the ind_pool...

            % log determinat:
            Det(k)=-subj_n/2*(log(1-p+bn1*p)+(bn1-1)*log(1-p)); % tested det
            TR(k)=-.5*(dgW*bn1-(sum(sum(W(ind_pool,ind_pool))))*(p/(1-p)/(1-p+bn1*p)));
            %toc (1/(1-p)).*(eye(bn1)-p/(1-p+bn1*p))
            %inv((1-p).*eye(bn1)+p)
            K_sets=[K_sets, ind_pool];

        end % each luster
        % then the edges between clusters...
        
        ind_btw=find(ismember(1:length(I),K_sets)==0);
        lik_btw=-0.5*length(ind_btw)*W(1,1);
        Lik(c+1)=sum(Det)+sum(TR)+lik_btw;
        
        Lik_norm=Lik-median(Lik);
        [lik_norm_max,lik_norm_maxid]=max(Lik-median(Lik))
        if lik_norm_max>100
           K(g)=K_label(lik_norm_maxid);
           
        else
            P=exp(Lik_norm)/sum(exp(Lik_norm))
            %P(end)=P(end)/10000;
            K(g)=K_label(mnrnd(1,P)==1); % label of cluster;
        end
         tabulate(K);
        if m>0, 
        para.K(m,:,g)=K;
        para.lik(m,g)=max(Lik);
        end
    end % for regions assigned g 
         % mcmc
         %tabulate(K)
%         if m>500 
%         para.K(m-500,g)=K;
%         para.lik(m-500,g)=lik_norm_max;
%         end
        
        
end


 
%% MCMC diagnostics
%para90=para; 
%para1=para;
% different starting points

%likelihood 
likh=(reshape(para.lik',numel(para.lik),1));
figure;plot(likh) 
%number of networks
 Kc1=[];
 for i=1:M
     for g=1:G
         Kc1=[Kc1,length(unique(para.K(i,:,g)))];
     end
 end
figure;         
 plot(Kc1)
% coverge
kmap=[[1:G]', squeeze(para.K(1,:,:)) squeeze(para.K(2,:,:)) squeeze(para.K(3,:,:)) squeeze(para.K(4,:,:))]; 
       % matrix to draw
 figure;
 imagesc(kmap);        % draw image and scale colormap to values range
 colorbar;          % show color scale


 kmap=[ones(G,1), squeeze(para.K(1,:,:)) squeeze(para.K(2,:,:) ) squeeze(para.K(3,:,:)) squeeze(para.K(4,:,:))]; 
       % matrix to draw
 figure;
 imagesc(kmap);        % draw image and scale colormap to values range
 colorbar;          % show color scale
 
 
  %% histogram
 unique(K)
ind_poolk=[]; % search index of the edges in the same cluster
            for i=1:length(I),

                w1=prod(single(ismember([I(i),J(i)],find(K==10))));
                 w2=prod(single(ismember([I(i),J(i)],find(K==24))));
                         w2=prod(single(ismember([I(i),J(i)],find(K==28))));
                if  w1%w1+w2+w3
                    ind_poolk=[ind_poolk,i]; % for example, 50 nodes will include 1225 edges in the communities
                end

            end
 
  pconcor_safe_in=1-pdist(X(ind_poolk,:),'correlation');      
   



 pconcor_safe_out=1-pdist(X(setdiff([1:435],ind_poolk),:),'correlation'); 
 
 figure; h1=histogram([pconcor_safe_out],40);
  hold on;
  h2=histogram([pconcor_safe_in]-.2,40);
  h1.Normalization = 'probability';
h1.BinWidth = 0.025;
h2.Normalization = 'probability';
h2.BinWidth = 0.025;