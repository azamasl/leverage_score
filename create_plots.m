function create_plots()

load('allresults');

for sn = 1:length(allresults.input.sampleNum)
    unif_err(sn) = allresults.unif_Data(sn).err;
    exac_err(sn) = allresults.exac_Data(sn).err;
    spec_err(sn) = allresults.spec_Data(sn).err;
    powe_err(sn) = allresults.powe_Data(sn).err;
    frob_err(sn) = allresults.frob_Data(sn).err;
    
    unif_time(sn) = allresults.unif_Data(sn).timings;
    exac_time(sn) = allresults.exac_Data(sn).timings; 
    spec_time(sn) = allresults.spec_Data(sn).timings;
    powe_time(sn) = allresults.powe_Data(sn).timings;   
    frob_time(sn) = allresults.frob_Data(sn).timings; 
end
%% Approximation error plots
fsize = 16;
ms = 10;
figure();
plot( allresults.input.sampleNum , unif_err/allresults.err, 'rd-','markers', ms);
hold on;
plot( allresults.input.sampleNum , exac_err/allresults.err, 'ks-','markers', ms);
hold on;
plot( allresults.input.sampleNum , frob_err/allresults.err, 'm^-','markers', ms);
hold on;
plot( allresults.input.sampleNum , spec_err/allresults.err, 'bo-','markers', ms);
hold on;
plot( allresults.input.sampleNum , powe_err/allresults.err, 'g*-','markers', ms);
hold off;

tit = 'Comparing approximation error';
legend({'Uniform','Exact','Frobenius', 'Spectral', 'Power method'}, 'FontSize', fsize)
title(tit,'FontSize',fsize)
ylabel('Error', 'FontSize', fsize)
xlabel('Number of sampling', 'FontSize', fsize)


%% Timing plots
figure();
semilogy( allresults.input.sampleNum , unif_time, 'rd-','markers', ms);
hold on;
semilogy( allresults.input.sampleNum , exac_time, 'ks-','markers', ms);
hold on;
semilogy( allresults.input.sampleNum , frob_time, 'm^-','markers', ms);
hold on;
semilogy( allresults.input.sampleNum , spec_time, 'bo-','markers', ms);
hold on;
semilogy( allresults.input.sampleNum , powe_time, 'g*-','markers', ms);
hold off;

legend({'Uniform','Exact','Frobenius', 'Spectral', 'Power method'}, 'FontSize', fsize)
title('Comparing the running time','FontSize',fsize)
ylabel('Time (s)', 'FontSize', fsize)
xlabel('Number of sampling', 'FontSize', fsize)
%close all
end