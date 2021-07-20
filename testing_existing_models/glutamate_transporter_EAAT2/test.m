[Glu_out_data, I_Glu_data] = csvimport('data/levy_etal_1998_fig2b.csv', 'columns', [1, 2] ,'noHeader', true);
Glu_out_data = Glu_out_data*1e-3; %convert from microM to mM

Glu_out = logspace(log10(min(Glu_out_data)),log10(max(Glu_out_data)),20);
kglu = 1.5e-2;
nglu = 2;
test_model = @(Glu) (Glu.^nglu./(kglu.^nglu+Glu.^nglu));
semilogx(Glu_out, test_model(Glu_out)); hold on
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, normalized');
set(gca,'FontSize',20);
semilogx(Glu_out_data, I_Glu_data, 'b^','MarkerSize',10,'LineWidth',2)


test2 = @(x) 0.0134.*(x.^nglu)./(kglu.^nglu+x.^nglu).*(-110-(13.35*log(1.8e5*x))).*exp(-0.0165*(-110-13.35*log(1.8e5*x)));
semilogx(Glu_out, -test2(Glu_out)./max(abs(test2(Glu_out))));

figure(2); semilogx(Glu_out, test2(Glu_out))