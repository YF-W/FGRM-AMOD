clear
load("Handwritten_5_258_6view");
 num = length(X); 
auc = zeros(1, num); 

tic;
for i = 1:num     
    Xs = X{i}; 
    gnd = out_label{i};
    s = FRGM_AMOD(Xs);
    [~, ~, ~, auc(i)] = perfcurve(gnd, s, 1);
end
elapsedTime = toc;

disp('Mean AUC:');
disp(mean(auc, 2)');
disp('Elapsed Time:');
disp(elapsedTime);

