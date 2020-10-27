function [ModelOutput] = BathMat_ShortTermModel(Chem,Modellength,Conc,TreatmentConc,Umean,SiteDepth,Dist2Shore,CageVolume)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

T = Modellength; % hours
t=T*3600;
k = 0.1;

%Plume Advction
longAdvec = Umean*t;
latAdvec = 4*(2*k*t)^0.5;
mixingDepth = min([SiteDepth/2, 10]);


if latAdvec < Dist2Shore*1000
    Vol = pi*(2*k*t)^0.5*longAdvec*mixingDepth;
elseif latAdvec > Dist2Shore*1000
    Vol = 0.5*pi*(Dist2Shore*1000)*longAdvec*mixingDepth;
end
%v1 = 0.5*pi*(Dist2Shore*1000)*longAdvec*mixingDepth;
%v2 = pi*(2*k*t)^0.5*longAdvec*mixingDepth;
%Vol = min([v1, v2]);
Area = Vol/mixingDepth;

consentMass  = Vol*Conc/1000000;
meanConc = CageVolume*TreatmentConc/(Area*mixingDepth);
peakConc = meanConc*(1/0.6);
treatmentVol = Area*mixingDepth*Conc/TreatmentConc;

noCagesTreated = treatmentVol/CageVolume;
areaExceedsEQS = 0.5*Area*0.000001;

ModelOutput = table();
ModelOutput.Chemical = {Chem};
ModelOutput.Time = T;                              %hrs
ModelOutput.distanceFromCagetoShore = Dist2Shore;  %m
ModelOutput.diffusionCoefficient = k;              %m^2/s
ModelOutput.Umean =Umean;                          %m/s
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

