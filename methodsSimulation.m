%% ode cooperation evaluation

n_genes = 100;  % Number of genes in the GRN

t_span = [0 20];                % Time interval for ODE simulation
k = 3;  % Number of expected connections per gene
wself = rand(1,n_genes) * 2;
x0 = rand(1, n_genes);

% Creating a base model
w3_ = rand(n_genes,n_genes,n_genes) * 2;            % The weights for each gene (range: 0–2)
L3 = randi(n_genes^2,n_genes,n_genes,n_genes)<=k;    % The structure connections between the genes
w3 = w3_.*L3;                   % Apply base weights to the base structure

ode1 = @(t, x, w3, w2) -x + coop(x,w3); % basic
ode2 = @(t, x, w3, w2) -x + w2'.*(x ./ (1 + x)) + coop(x,w3); % with self activation
ode3 = @(t, x, w3, w2) -x + w2'.*(x ./ (1 + x)); % only self activation

[t1, x1] = ode45(@ (t, x) ode1(t,x,w3,wself), t_span, x0); % Solve ODE
[t2, x2] = ode45(@ (t, x) ode2(t,x,w3,wself), t_span, x0); % Solve ODE
[t3, x3] = ode45(@ (t, x) ode3(t,x,w3,wself), t_span, x0); % Solve ODE

   

%plotting
figure("Name","basic"); 
plot(t1,x1,'-o')
title("basic cooperation")
xlabel("step")
ylabel("expression levels")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

figure("Name","with self-activation"); 
plot(t2,x2,'-o')
title("cooperation with self-activation")
xlabel("step")
ylabel("expression levels")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

figure("Name","only self-activation"); 
plot(t3,x3,'-o')
title("only self-activation")
xlabel("step")
ylabel("expression levels")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);


%% p effect model 1

n_genes = 1000;  % Number of genes in the GRN
n_cells = 100;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';

t_span = [0 20];                % Time interval for ODE simulation
wself = rand(1,n_genes) * 2;

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1 = randi(n_genes,n_genes)<=k; % The structure connections between the genes
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

data0 = zeros(n_genes,n_cells);  % final simulated data for p=0
p = 0;
    % creating the data for p=0
for i = 1:n_cells
    L2 = rand(n_genes)<=p;          % Random mask
    L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
    L2 = logical(L2.*L);            % Apply the mask on the base structure
    
    w = w_.*L;                      % Apply base weights to the base structure
    w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
    
    x0 = rand(1, n_genes);
    offgenes = randi(n_genes,1,5);
    x0(offgenes) = 0;
    w(offgenes,:) = 0;
    w(:,offgenes) = 0;

    [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
    data0(:,i) = x(end,:); % Save final state
end
data0(data0<10^-7)=0;

data1 = zeros(n_genes,n_cells);  % final simulated data for p=1
p = 1;
    % creating the data for p=1
for i = 1:n_cells
    L2 = rand(n_genes)<=p;          % Random mask
    L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
    L2 = logical(L2.*L);            % Apply the mask on the base structure
    
    w = w_.*L;                      % Apply base weights to the base structure
    w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
    
    x0 = rand(1, n_genes);
    offgenes = randi(n_genes,1,5);
    x0(offgenes) = 0;
    w(offgenes,:) = 0;
    w(:,offgenes) = 0;

    [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
    data1(:,i) = x(end,:); % Save final state
end
data1(data1<10^-7)=0;

%plotting
figure("Name","p effect model 1"); % open new figure named the type of data
hold all;

%plooting p=0
subplot(1,2,1);
imagesc(data0)
title("p = 0")
xlabel("cells")
ylabel("genes")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plooting p=1
subplot(1,2,2);
imagesc(data1)
title("p = 1")
xlabel("cells")
ylabel("genes")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);


%% comparison with p change model 1

n_genes = 1000;  % Number of genes in the GRN
n_cells = 100;
p = 0.1:0.1:1;
alpha = 10.^-4;
k = 3;  % Number of expected connections per gene

ode1 = @(t, x, w, w2) -x + w2'.*(x ./ (1 + x)) + sum(w.*(x ./ (1 + x)),1)';


% bootstrap characteristics
numsub = 30; % number of times
persub = 0.75; % precent from the data

PD = zeros(length(p),2);
precision = zeros(length(p),2);
results = zeros(length(p),2);


t_span = [0 20];                % Time interval for ODE simulation
wself = rand(1,n_genes) * 2;

% Creating a base model
w_ = rand(n_genes) * 2;         % The weights for each gene (range: 0–2)
L1 = randi(n_genes,n_genes)<=k; % The structure connections between the genes
L = triu(L1,1)+triu(L1,1)';     % Make the matrix symmetric

data = zeros(n_genes,n_cells);  % Matrix to store final simulated data

for j = 1:length(p)
    % creating the data
    for i = 1:n_cells
        L2 = rand(n_genes)<=p(j);          % Random mask
        L2 = triu(L2,1)+triu(L2,1)';    % Symmetrize
        L2 = logical(L2.*L);            % Apply the mask on the base structure
        
        w = w_.*L;                      % Apply base weights to the base structure
        w(L2) = rand(sum(L2,"all"),1) * 2; % Randomize weights according to mask
        
        x0 = rand(1, n_genes);
        offgenes = randi(n_genes,1,5);
        x0(offgenes) = 0;
        w(offgenes,:) = 0;
        w(:,offgenes) = 0;
    
        [t, x] = ode45(@ (t, x) ode1(t,x,w,wself), t_span, x0); % Solve ODE
        data(:,i) = x(end,:); % Save final state
    end
    data(data<10^-7)=0;

    BSgcn = zeros(1,numsub); % the vector of results
    confo = zeros(2,numsub);

    parfor f=1:numsub
        index = randperm(n_cells); % randomly picking the cells.
        data2 = data(:,index(1:floor(n_cells*persub))); % lessen the data.
        [BSgcn(f),A] = GCN(data2,alpha); % perform the function
        compatibility = (A-1)*2+L;
        TP = sum(compatibility==1,'all');
        FN  = sum(compatibility==-1,'all');
        FP = sum(compatibility==0,'all');
        
        confo(:,f) = [TP./(TP + FN),TP./(TP + FP)];
    end

    results(j,:) = [mean(BSgcn) std(BSgcn)] ;

    PD(j,:) = [mean(confo(1,:)) std(confo(1,:))];
    precision(j,:) = [mean(confo(2,:)) std(confo(2,:))];
  
end

%plotting
figure("Name",n_cells+"X"+n_genes); % open new figure named the type of data
hold all;

%plotting the gcn
subplot(2,1,1);
errorbar(p,results(:,1),results(:,2),'LineWidth',2)
title("MND")
xlabel("p")
ylabel('MND')
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plooting the recall
subplot(2,2,3);
errorbar(p,PD(:,1),PD(:,2),'LineWidth',2)
title("Recall")
xlabel("p")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);

%plotting the precision
subplot(2,2,4);
errorbar(p,precision(:,1),precision(:,2),'LineWidth',2)
title("Precision")
xlabel("p ")
set(gcf, 'Color', 'w');
box on;
set(gca, 'FontSize', 16);