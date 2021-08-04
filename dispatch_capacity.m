%%%demand side
zerocost = 7504.1;%0cost generators including nuclear hydro and bio:ercot = 7504.1 pjm = 40667
load = load1;
netloadnw = load(:,1) - zerocost; %netload without wind
wexp = windall;% wind forecast expecation
wsigma = windallsigma;% wind forecast standard deviation
winmax = wexp + 3 * wsigma;
netload = netloadnw - wexp;

h = size(load);
h = h(1);
n = size(mc);
n = n(1);

%%%solve for the DA dispatch_static dispatch
upbound = eye(n);
loadbalance = -1*ones(1,n);
A = [upbound;loadbalance];
netload = -1*netload;
netloadnw = -1*netloadnw;
SIdx = zeros(h);% index of the marginal generator for hr i
solution = zeros(n,h);
solutionnw = zeros(n,h);
pda = zeros(h,1);%day-ahead price
mneda = zeros(h,1);%day-ahead marginal nox emission
mseda = zeros(h,1);%day-ahead marginal nox emission
mceda = zeros(h,1);%day-ahead marginal nox emission
mcheda = zeros(h,1);%day-ahead marginal nox emission
lb = zeros(n,1);
for i = 1:h
    [x,fval,exitflag,output,lambda] = linprog(mc,loadbalance,netload(i),[],[],lb,capacity);% dispatch with wind
    [xnw,fvalnw,exitflagnw,outputnw,lambdanw] = linprog(mc,loadbalance,netloadnw(i),[],[],lb,capacity);%dispatch wo wind
    solution(:,i) = x;
    solutionnw(:,i) = xnw;
    index = find(x>0.1);
    SIdx(i) = length(index);
    resid(:,i) = ramp;
    j = max(SIdx(i),1);
    resid(j,i) = min(capacity(j) - solution(j,i), ramp(j));
    pda(i) = lambda.ineqlin;
    mneda(i) = ner(j);
    mseda(i) = ser(j);
    mceda(i) = cer(j);
    mcheda(i) = cher(j);
end


%ERCOT / PJM
%%%calculate the expectation of emissoin abatement from integrating wind  
%%%energy
tcexp = zeros(h,1);%expectation of real-time dispatch total cost
neexp = zeros(h,1);%expectation of emission from wind shortfalls
seexp = zeros(h,1);
ceexp = zeros(h,1);
cheexp = zeros(h,1);

mcrtexp = zeros(h,1);%real-time expected marginal cost
mnert = zeros(h,1);%marginal expected emission from wind shortfalls
msert = zeros(h,1);
mcert = zeros(h,1);
mchert = zeros(h,1);
rtgenexp = zeros(h,1);%expected total generation in RT process

RTload = wexp; % the load need to be balanced in the real-time market, in the tradtional situation, this is the wind forecast value.


SIdxrt = zeros(h);% index of the marginal generator for hr i in real-time market
for i = 1:h
    j = max(SIdx(i),1);
    r = RTload(i);
  
    z = normcdf(winmax(i),wexp(i),wsigma(i)) - normcdf(0,wexp(i),wsigma(i));% the parameter for calculaing the truncated normal CDF and PDF
    while r - resid(j,i)> 0 %if the demand needed to be balanced is larger than the available capaicity of generator i in hour j in the RT market, the step width is the available capacity of generator i
        prb = (normcdf(r,wexp(i),wsigma(i))-normcdf(r - resid(j,i),wexp(i),wsigma(i)))/z;
        fgenexp = @(w)(r-w) .* normpdf(w,wexp(i),wsigma(i))/z; 
        genexp = quad(fgenexp,r-resid(j,i),r);
        rtgenexp(i) = rtgenexp(i) + genexp; 
        tcexp(i) = tcexp(i) + mc(j) * genexp;
        neexp(i) = neexp(i) + ner(j) * genexp; 
        seexp(i) = seexp(i) + ser(j) * genexp;
        ceexp(i) = ceexp(i) + cer(j) * genexp;
        cheexp(i) = cheexp(i) + cher(j) * genexp;
        mcrtexp(i) = mcrtexp(i) + mc(j) * prb;
        mnert(i) =  mnert(i)+ ner(j) * prb;
        msert(i) =  msert(i)+ ser(j) * prb;
        mcert(i) =  mcert(i)+ cer(j) * prb;
        mchert(i) = mchert(i)+ cher(j) * prb;
        r = r - resid(j,i);
        j = j+1;
    end

% calculate the effects of last generator if the demand need to be balanced
% is smaller than the available capaicity of generator i in hour j in the
% RT market, the step width is the residual demand to be balanced. 
  
    prb = (normcdf(r,wexp(i),wsigma(i))-normcdf(0,wexp(i),wsigma(i)))/(normcdf(winmax(i),wexp(i),wsigma(i)) - normcdf(0,wexp(i),wsigma(i)));
    fgenexp = @(w)(r-w) .* normpdf(w,wexp(i),wsigma(i))/z; 
    genexp = quad(fgenexp,0,r);
    rtgenexp(i) = rtgenexp(i) + genexp; 
    rtgenexp(i) = rtgenexp(i) + genexp; 
        tcexp(i) = tcexp(i) + mc(j) * genexp;
        neexp(i) = neexp(i) + ner(j) * genexp; 
        seexp(i) = seexp(i) + ser(j) * genexp;
        ceexp(i) = ceexp(i) + cer(j) * genexp;
        cheexp(i) = cheexp(i) + cher(j) * genexp;
        mcrtexp(i) = mcrtexp(i) + mc(j) * prb;
        mnert(i) =  mnert(i)+ ner(j) * prb;
        msert(i) =  msert(i)+ ser(j) * prb;
        mcert(i) =  mcert(i)+ cer(j) * prb;
        mchert(i) = mchert(i)+ cher(j) * prb;
       r = 0;  
    SIdxrt(i)=j;
end



%%%calculate emission abatement and expected wind generation
noxnw = ner'*solutionnw;
so2nw = ser'*solutionnw;
co2nw = cer'*solutionnw;
ch4nw = cher'*solutionnw;
ttgcostnw = mc'* solutionnw ; 

noxw = ner'*solution + neexp';
so2w = ser'*solution + seexp';
co2w = cer'*solution + ceexp';
ch4w = cher'*solution + cheexp';
ttgcost = mc'* solution + tcexp'; 

noxabt = noxnw - noxw;
so2abt = so2nw - so2w;
co2abt = co2nw - co2w;
ch4abt = ch4nw - ch4w;
costabt = ttgcostnw-ttgcost;

Capa = [0.45359237*co2abt'./wexp, 0.45359237*so2abt'./wexp, 0.45359237*noxabt'./wexp, 0.45359237*ch4abt'./wexp, costabt'./wexp];

Capa_all = [max(Capa);min(Capa);mean(Capa, 1);std(Capa, 1);];

Endowment = Capa./repmat(max(Capa),[size(Capa,1),1]);
Ave_endow = mean(Endowment, 2);
com_advan = Endowment./repmat(Ave_endow, [1, size(Endowment, 2)]);

%ERCOT two region / PJM two region 
%%%calculate the expectation of emissoin abatement from integrating wind  
%%%energy
wexp = windnorth;% wind forecast expecation
wsigma = windnorthsigma;% wind forecast standard deviation
rwexp = windsouth;% wind forecast expecation
winmax = wexp + 3 * wsigma;
netloadnw = load(:,1) - zerocost; %netload without wind
netload = netloadnw - wexp;

%%%solve for the DA dispatch_static dispatch
netload = -1*netload;
netloadnw = -1*netloadnw;
SIdx = zeros(h);% index of the marginal generator for hr i
solution = zeros(n,h);
solutionnw = zeros(n,h);
pda = zeros(h,1);%day-ahead price
mneda = zeros(h,1);%day-ahead marginal nox emission
mseda = zeros(h,1);%day-ahead marginal nox emission
mceda = zeros(h,1);%day-ahead marginal nox emission
mcheda = zeros(h,1);%day-ahead marginal nox emission
lb = zeros(n,1);
for i = 1:h
    [x,fval,exitflag,output,lambda] = linprog(mc,loadbalance,netload(i),[],[],lb,capacity);% dispatch with wind
    solution(:,i) = x;
    index = find(x>0.1);
    SIdx(i) = length(index);
    resid(:,i) = ramp;
    j = max(SIdx(i),1);
    resid(j,i) = min(capacity(j) - solution(j,i), ramp(j));
    pda(i) = lambda.ineqlin;
    mneda(i) = ner(j);
    mseda(i) = ser(j);
    mceda(i) = cer(j);
    mcheda(i) = cher(j);
end


tcexp = zeros(h,1);%expectation of real-time dispatch total cost
neexp = zeros(h,1);%expectation of emission from wind shortfalls
seexp = zeros(h,1);
ceexp = zeros(h,1);
cheexp = zeros(h,1);

mcrtexp = zeros(h,1);%real-time expected marginal cost
mnert = zeros(h,1);%marginal expected emission from wind shortfalls
msert = zeros(h,1);
mcert = zeros(h,1);
mchert = zeros(h,1);
rtgenexp = zeros(h,1);%expected total generation in RT process

RTload = wexp; % the load need to be balanced in the real-time market, in the tradtional situation, this is the wind forecast value.


SIdxrt = zeros(h);% index of the marginal generator for hr i in real-time market
for i = 1:h
    j = max(SIdx(i),1);
    r = RTload(i);
  
    z = normcdf(winmax(i),wexp(i),wsigma(i)) - normcdf(0,wexp(i),wsigma(i));% the parameter for calculaing the truncated normal CDF and PDF
    while r - resid(j,i)> 0 %if the demand needed to be balanced is larger than the available capaicity of generator i in hour j in the RT market, the step width is the available capacity of generator i
        prb = (normcdf(r,wexp(i),wsigma(i))-normcdf(r - resid(j,i),wexp(i),wsigma(i)))/z;
        fgenexp = @(w)(r-w) .* normpdf(w,wexp(i),wsigma(i))/z; 
        genexp = quad(fgenexp,r-resid(j,i),r);
        rtgenexp(i) = rtgenexp(i) + genexp; 
        tcexp(i) = tcexp(i) + mc(j) * genexp;
        neexp(i) = neexp(i) + ner(j) * genexp; 
        seexp(i) = seexp(i) + ser(j) * genexp;
        ceexp(i) = ceexp(i) + cer(j) * genexp;
        cheexp(i) = cheexp(i) + cher(j) * genexp;
        mcrtexp(i) = mcrtexp(i) + mc(j) * prb;
        mnert(i) =  mnert(i)+ ner(j) * prb;
        msert(i) =  msert(i)+ ser(j) * prb;
        mcert(i) =  mcert(i)+ cer(j) * prb;
        mchert(i) = mchert(i)+ cher(j) * prb;
        r = r - resid(j,i);
        j = j+1;
    end

% calculate the effects of last generator if the demand need to be balanced
% is smaller than the available capaicity of generator i in hour j in the
% RT market, the step width is the residual demand to be balanced. 
  
    prb = (normcdf(r,wexp(i),wsigma(i))-normcdf(0,wexp(i),wsigma(i)))/(normcdf(winmax(i),wexp(i),wsigma(i)) - normcdf(0,wexp(i),wsigma(i)));
    fgenexp = @(w)(r-w) .* normpdf(w,wexp(i),wsigma(i))/z; 
    genexp = quad(fgenexp,0,r);
    rtgenexp(i) = rtgenexp(i) + genexp; 
    rtgenexp(i) = rtgenexp(i) + genexp; 
        tcexp(i) = tcexp(i) + mc(j) * genexp;
        neexp(i) = neexp(i) + ner(j) * genexp; 
        seexp(i) = seexp(i) + ser(j) * genexp;
        ceexp(i) = ceexp(i) + cer(j) * genexp;
        cheexp(i) = cheexp(i) + cher(j) * genexp;
        mcrtexp(i) = mcrtexp(i) + mc(j) * prb;
        mnert(i) =  mnert(i)+ ner(j) * prb;
        msert(i) =  msert(i)+ ser(j) * prb;
        mcert(i) =  mcert(i)+ cer(j) * prb;
        mchert(i) = mchert(i)+ cher(j) * prb;
       r = 0;  
    SIdxrt(i)=j;
end

rnoxw = ner'*solution + neexp';
rso2w = ser'*solution + seexp';
rco2w = cer'*solution + ceexp';
rch4w = cher'*solution + cheexp';
rttgcost = mc'* solution + tcexp'; 

noxabt = rnoxw - noxw;
so2abt = rso2w - so2w;
co2abt = rco2w - co2w;
ch4abt = rch4w - ch4w;
costabt = rttgcost - ttgcost ; 

Capa = [0.45359237*co2abt'./rwexp, 0.45359237*so2abt'./rwexp, 0.45359237*noxabt'./rwexp, 0.45359237*ch4abt'./rwexp, costabt'./rwexp];

Capa_all = [max(Capa);min(Capa);mean(Capa, 1);std(Capa, 1);];

Endowment = Capa./repmat(max(Capa),[size(Capa,1),1]);
Ave_endow = mean(Endowment, 2);
com_advan = Endowment./repmat(Ave_endow, [1, size(Endowment, 2)]);

%¾ÛÀà
j = 1;
b = zeros(30,1);
for a = 10:40
[Idx,Ctrs,sumD,D] = kmeans(comadvan, 25,'dist','sqEuclidean','Replicates',3);
x = 0;
for i = 1:length(comadvan)
    if D(i,Idx(i)) > 0.05
        x = x + 1;
    end
end 
b(j,1) = x/length(comadvan);
j = j + 1;
end
  
    
