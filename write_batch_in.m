function [output] = write_batch_in(O2, T, P, phi, t)
% This function creates the input for the batch

N2 = 1 - O2;

fid = fopen('input.dic', 'w');
fprintf(fid, 'Dictionary BatchReactor \n { \n');
fprintf(fid, '@KineticsPreProcessor	   kinetic-mechanism; \n');
fprintf(fid, '@Type					   NonIsothermal-ConstantVolume; \n');
fprintf(fid, '@InitialStatus          initial-mixture; \n');
% fprintf(fid, '@SensitivityAnalysis     sensitivity; \n');
fprintf(fid, '@OnTheFlyROPA      ropa; \n');
fprintf(fid, '@OnTheFlyPostProcessing post-processing; \n');
fprintf(fid, '@EndTime				%d s; \n } \n', 0.01);

fprintf(fid, 'Dictionary kinetic-mechanism \n { \n');
fprintf(fid, '@Kinetics           ../aramco/chemAramco.cki; \n');
fprintf(fid, '@Thermodynamics		../aramco/thermAramco.ckt; \n');
fprintf(fid, '@Output             kinetics; \n } \n');

fprintf(fid, 'Dictionary initial-mixture \n { \n');
fprintf(fid, '@Temperature   	%d   K; \n', T);
fprintf(fid, '@Pressure      	%d Pa; \n', P);
fprintf(fid, '@OxidizerMoleFractions        O2 %d \n', O2);
fprintf(fid, '                              N2 %d; \n', N2);
fprintf(fid, '@FuelMoleFractions   CH4 1; \n');
fprintf(fid, '@EquivalenceRatio     %d; \n } \n', phi);

fprintf(fid, 'Dictionary sensitivity \n { \n');
fprintf(fid, '     @Type   arrhenius-parameters; \n');
fprintf(fid, '    @DenseSolver    Eigen; \n');
fprintf(fid, '@DenseFullPivoting     false; \n');
fprintf(fid, '     @SubSteps     5; \n');
fprintf(fid, '     @Species     OH O2 CH4 CH3; \n } \n');

fprintf(fid, 'Dictionary    ropa \n { \n');
fprintf(fid, '@Species     CH4 OH O2 H2O2 HO2; \n');
fprintf(fid, '@ReferenceSpecies CH4; \n');
fprintf(fid, '@Threshold 0.01; \n');
fprintf(fid, '@Times     %d s; \n }  \n', t);

fprintf(fid, 'Dictionary    post-processing \n { \n');
fprintf(fid, '@ReactionRates     129; \n } \t');

output = true;


end

