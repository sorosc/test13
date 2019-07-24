r_VBISX = price2ret(VBISX);
r_VBLTX = price2ret(VBLTX);
r_VEIEX = price2ret(VEIEX);
r_VEURX = price2ret(VEURX);
r_VFINX = price2ret(VFINX);
returns = [r_VFINX r_VEURX r_VEIEX r_VBLTX r_VBISX];
N = 5;
%3a
SIGMA = cov(returns);
w0 = ones(N,1)/N;
Aeq = ones(1,N);
beq = 1;
%wa = fmincon(@(w) w'*SIGMA*w,w0,[],[],Aeq,beq,-ones(N,1),ones(N,1));
wa = fmincon(@(w) w'*SIGMA*w,w0,[],[],Aeq,beq);
means = mean(returns,1);
Era = means*wa;
Sda = sqrt(wa'*SIGMA*wa); 
display(['Part a) expected return = ' num2str(100*Era) '%']);
display(['                     SD = ' num2str(100*Sda) '%']);
figure;
bar(wa,'k');
xlabel('Asset'); ylabel('Weight');
%3b
wb = fmincon(@(w) w'*SIGMA*w,w0,[],[],Aeq,beq,zeros(N,1),ones(N,1));
Erb = means*wb;
Sdb = sqrt(wb'*SIGMA*wb); 
display(['Part b) expected return = ' num2str(100*Erb) '%']);
display(['                     SD = ' num2str(100*Sdb) '%']);
figure;
bar(wb,'k');
xlabel('Asset'); ylabel('Weight');
%3c
display('Part c)');
% p = Portfolio('assetmean', means*100, 'assetcovar', SIGMA*100, ...
% 'lowerbudget', 1, 'upperbudget', 1, 'lowerbound', -1);
% p = estimateAssetMoments(p,returns*100,'missingdata',false);
% Plot Efficient Frontier
A = ones(2,N);
A(2,:) = means;
M = 100;
step = (max(means)-Era)/M;
rmesh = Era:step:max(means);
sdmesh = nan(1,M+1);
for i=1:M+1
    b = [1;rmesh(i)];
    w = inv(SIGMA)*A'*inv(A*inv(SIGMA)*A')*b;
    sdmesh(i) = sqrt(w'*SIGMA*w);
end
figure;
plot(sdmesh*100,rmesh*100,'k');
% Plot the lower half
rmesh2 = Era-step*(M/5+1):step:Era;
sdmesh2 = nan(1,M/5);
for i=1:(M/5+2)
    b = [1;rmesh2(i)];
    w = inv(SIGMA)*A'*inv(A*inv(SIGMA)*A')*b;
    sdmesh2(i) = w'*SIGMA*w;
end
xlabel('Portfolio SD (%)'); ylabel('Portfolio Monthly Return (%)');
% M = 10;
% pwgt = estimateFrontier(p,M);
% pnames = cell(1,M);
% for i = 1:M
%     pnames{i} = sprintf('Port%d',i);
% end
AssetList = {'S&P 500','Euro stock','Emerging market','L/T bond','S/T bond','Pacific stock'};
% p = Portfolio('AssetList',AssetList);
% Blotter = dataset([{pwgt},pnames],'obsnames',p.AssetList);
% disp(Blotter);
%3d
rf = 0.0004167;
slopes = (rmesh-rf)./sdmesh;
[slope,maxid] = max(slopes);
bestr = rmesh(maxid);
bestsd = sdmesh(maxid);
SD_assets = zeros(5,1);
for i =1:5
    SD_assets(i,1) = sqrt(SIGMA(i,i));
end
figure;
plot(sdmesh*100,rmesh*100,'k');
hold on;
plot(sdmesh2*100,rmesh2*100,'r');
xlabel('Portfolio SD (%)'); ylabel('Portfolio Monthly Return (%)');
hold on;
plot([0 bestsd]*100,[rf bestr]*100,'k--');
title('Part d)');
text(100*bestsd,100*bestr,['\leftarrow Tangency Portfolio (Return = ' ...
num2str(bestr*100) '% SD = ' num2str(bestsd*100) '%)'],'HorizontalAlignment','left');
display(['Part d) tangency return = ' num2str(100*bestr) '%']);
display([' SD = ' num2str(100*bestsd) '%']);
hold on;
scatter(SD_assets*100,means*100);
text(SD_assets(1,1)*100+0.1,means(1,1)*100,['S&P 500'],'HorizontalAlignment','left');
text(SD_assets(2,1)*100-0.1,means(1,2)*100,['Euro Stock'],'HorizontalAlignment','right');
text(SD_assets(3,1)*100-0.1,means(1,3)*100,['Emerging Market'],'HorizontalAlignment','right');
text(SD_assets(4,1)*100+0.1,means(1,4)*100,['Long-Term Bond'],'HorizontalAlignment','left');
text(SD_assets(5,1)*100+0.1,means(1,5)*100,['Short-Term Bond'],'HorizontalAlignment','left');
text(100*Sda,100*Era+0.01,['\leftarrow Global Minimum Variance Portfolio'],'HorizontalAlignment','left');
text(0.04,0.07,['\leftarrow Capital Allocation Line (T-bills + Tangency Portfolio)'],'HorizontalAlignment','left');
text(0,rf*100,['\leftarrow Risk-Free Rate'],'HorizontalAlignment','left');
% Bar chart of weights
b = [1;bestr];
bestw = inv(SIGMA)*A'*inv(A*inv(SIGMA)*A')*b;
figure;
bar(bestw,'k');
title('Part d) Tangency Portfolio');
xlabel('Asset'); ylabel('Weight');
% Sharpe ratio for each asset
display('Sharpe Ratio of each asset');
for i=1:N
    display([AssetList{i} '=' num2str((means(i)-rf)/sqrt(SIGMA(i,i)))]);
end
display('S/T bond has the highest Sharpe ratio.');
display(['Sharpe Ratio of tangency port = ' num2str(slope)]);
%3e
figure;
plot(sdmesh*100,rmesh*100,'k');
hold on;
sdmesh = nan(1,M+1);
A(2,:) = -means;
for i=1:M+1
    b = [1;-rmesh(i)];
    w = fmincon(@(w) w'*SIGMA*w,w0,A,b,Aeq,beq,zeros(N,1),ones(N,1));
    sdmesh(i) = sqrt(w'*SIGMA*w);
end
plot(sdmesh*100,rmesh*100,'--');
title('Part e)');
legend(' Short-sell','No short-sell');
xlabel('Portfolio SD (%)'); ylabel('Portfolio Monthly Return (%)');
%3e_rev
figure;
plot(sdmesh*100,rmesh*100,'k');
hold on;
sdmesh = nan(1,M+1);
A(2,:) = -means;
Aeq2 = [A;Aeq];
for i=1:M+1
    b = [1;-rmesh(i)];
    Beq2 = [b;beq];
    w = fmincon(@(w) w'*SIGMA*w,w0,[],[],Aeq2,Beq2,zeros(N,1),ones(N,1));
    sdmesh(i) = sqrt(w'*SIGMA*w);
end
plot(sdmesh*100,rmesh*100,'--');
title('Part e)');
legend(' Short-sell','No short-sell');
xlabel('Portfolio SD (%)'); ylabel('Portfolio Monthly Return (%)');