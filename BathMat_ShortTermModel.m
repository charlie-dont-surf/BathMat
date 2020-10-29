function [ModelOutput] = BathMat_ShortTermModel(Chem,Input,ChemEQS,CageVolume)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Input variable units
%Conc = ng/l
%TreatmentConc = ng/l
%Umean = m/s
%SiteDepth = m
%Dist2Shore = km
%CageVolume = m^3
%% 

T = ChemEQS.Time; % hours
t=T*3600;
k = 0.1;

%Plume Advction
longAdvec = Input.Umean*t;
latAdvec = 4*(2*k*t)^0.5;
mixingDepth = min([Input.SiteDepth/2, 10]);


if latAdvec < Input.Dist2Shore*1000
    Vol = pi*(2*k*t)^0.5*longAdvec*mixingDepth;
elseif latAdvec > Input.Dist2Shore*1000
    Vol = 0.5*pi*(Input.Dist2Shore*1000)*longAdvec*mixingDepth;
end
%v1 = 0.5*pi*(Dist2Shore*1000)*longAdvec*mixingDepth;
%v2 = pi*(2*k*t)^0.5*longAdvec*mixingDepth;
%Vol = min([v1, v2]);
Area = Vol/mixingDepth;
Vol_lts = (Area*mixingDepth)*1000; %convert volume to lts
CageVolume_lts = Input.CageVolume*1000; %convert volume to lts

consentMass  = Vol*ChemEQS.EQSconc/1000000; %g
meanConc = CageVolume_lts*ChemEQS.TreatmentConc/Vol_lts;
peakConc = meanConc*(1/0.6);
treatmentVol = (Area*mixingDepth)*ChemEQS.EQSconc/ChemEQS.TreatmentConc;

noCagesTreated = treatmentVol/CageVolume;
areaExceedsEQS = 0.5*Area*0.000001;

ModelOutput = table();
ModelOutput.Chemical = {Chem};
ModelOutput.Time = T;                              %hrs
ModelOutput.distanceFromCagetoShore = Input.Dist2Shore;  %m
ModelOutput.diffusionCoefficient = k;              %m^2/s
ModelOutput.Umean = Input.Umean;                          %m/s
ModelOutput.plumeLength = longAdvec;               %m
ModelOutput.plumeWidth = latAdvec;                 %m
ModelOutput.time = T;                              %hrs
ModelOutput.mixingZoneArea = Area;                 %m^2
ModelOutput.treatmentVol = treatmentVol;           %m^3
ModelOutput.noCagesTreated = noCagesTreated;       %
ModelOutput.meanConc = meanConc;                   %ng/l
ModelOutput.consentMass = consentMass;             %g
ModelOutput.peakConc = peakConc;                   %ng/l
ModelOutput.areaExceedsEQS = areaExceedsEQS;       %km^2


end

