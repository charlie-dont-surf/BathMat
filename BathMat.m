%% BathMat
%Tool for modelling bath treatments

%% Input parameters
SiteName = 'Chalmers Hope';
SiteCenter = [328735, 1001311];
NoCages = 12;
CageCirc = 120;

SiteDepth = 30;
Dist2Shore = 0.49;

TreatmentDepth = 10; %either cone depth or flat depth %3.4m for flat or 10m for cone
cageShape = 'cone'; %flat or cone or wellboat
wellboatCapacity = 1000; %m^3

%Flow Parameters
Umean = 0.122;
%Upara = 0.004;
%Vpara = 0.012;
%Uamp = 0.183;
%Vamp = 0.094;


%Initail concentrations
CageArea = pi*(CageCirc/(2*pi))^2;
switch cageShape
    case {'flat'}
        CageVolume = CageArea * TreatmentDepth;
    case{'cone'}
        CageVolume = CageArea * TreatmentDepth/3;
    case{'wellboat'}
        CageVolume = wellboatCapacity;
end  
CageVolumeWellBoat = wellboatCapacity;


%% Short-term Model 1
%Azamethiphos
AzaEQS.Time = 3;                %hrs
AzaEQS.EQSconc = 250;           %ng/l
AzaEQS.TreatmentConc = 100000; %ng/l
%Cypermethrin
CypEQS.Time = 3;                %hrs
CypEQS.EQSconc = 0.06;          %ng/l
CypEQS.TreatmentConc = 5000;    %ng/l
%Deltamethrin
DelEQS.Time = 3;                %hrs
DelEQS.EQSconc = 9;             %ng/l
DelEQS.TreatmentConc = 2000;    %ng/l

ModelOutput1 = BathMat_ShortTermModel('AZA',AzaEQS.Time, AzaEQS.EQSconc, AzaEQS.TreatmentConc, Umean,SiteDepth,Dist2Shore,CageVolume);
ModelOutput2 = BathMat_ShortTermModel('AZA-WellBoat',AzaEQS.Time, AzaEQS.EQSconc, AzaEQS.TreatmentConc, Umean,SiteDepth,Dist2Shore,CageVolumeWellBoat);
ModelOutput3 = BathMat_ShortTermModel('CYP',CypEQS.Time, CypEQS.EQSconc, CypEQS.TreatmentConc, Umean,SiteDepth,Dist2Shore,CageVolume);
ModelOutput4 = BathMat_ShortTermModel('CYP-WellBoat',CypEQS.Time, CypEQS.EQSconc, CypEQS.TreatmentConc, Umean,SiteDepth,Dist2Shore,CageVolumeWellBoat);
ModelOutput5 = BathMat_ShortTermModel('DEL',DelEQS.Time, DelEQS.EQSconc, DelEQS.TreatmentConc, Umean,SiteDepth,Dist2Shore,CageVolume);
ModelOutput6 = BathMat_ShortTermModel('DEL-WellBoat',DelEQS.Time, DelEQS.EQSconc, DelEQS.TreatmentConc, Umean,SiteDepth,Dist2Shore,CageVolumeWellBoat);

BathMat1 = [ModelOutput1;ModelOutput2;ModelOutput3;ModelOutput4;ModelOutput5;ModelOutput6];
head(BathMat1)

%Calculate treatable volume from 24hr Aza run
AzaTreatmentMass = 168.5;
AzatreatmentVol = AzaTreatmentMass*10;





%%
% T = 3; % hours
% t=T*3600;
% k = 0.1;
% 
% %Plume Advction
% longAdvec = Umean*t;
% latAdvec = 4*(2*k*t)^0.5;
% mixingDepth = min([SiteDepth/2, 10]);
% 
% v1 = 0.5*pi*(Dist2Shore*1000)*longAdvec*mixingDepth;
% v2 = pi*(2*k*t)^0.5*longAdvec*mixingDepth;
% Vol = min([v1, v2]);
% Area = Vol/mixingDepth;
% 
% consentMass  = Vol*EQS.conc/1000000;
% meanConc = CageVolume*EQS.TreatmentConc/(Area*mixingDepth);
% peakConc = meanConc*(1/0.6);
% treatmentVol = Area*mixingDepth*EQS.conc/EQS.TreatmentConc;
% 
% noCagesTreated = treatmentVol/CageVolume;
% areaExceedsEQS = 0.5*Area*0.000001;
% 
% %Populate table
% BathMat1 =table();
% BathMat1.distanceFromCagetoShore = Dist2Shore;  %m
% BathMat1.diffusionCoefficient = k;              %m^2/s
% BathMat1.Umean =Umean;                          %m/s
% BathMat1.plumeLength = longAdvec;               %m
% BathMat1.plumeWidth = latAdvec;                 %m
% BathMat1.time = T;                              %hrs
% BathMat1.mixingZoneArea = Area;                 %m^2
% BathMat1.treatmentVol = treatmentVol;           %m^3
% BathMat1.noCagesTreated = noCagesTreated;       %
% BathMat1.meanConc = meanConc;                   %ng/l
% BathMat1.consentMass = consentMass;             %g
% BathMat1.peakConc = peakConc;                   %ng/l
% BathMat1.areaExceedsEQS = areaExceedsEQS;       %km^2
%head(BathMat1)

%% Long-term Model (Gillibrand and Turrell 1999)
ModelType = 'OpenWater';
mixingDepth = min([SiteDepth/2, 10]);
D = 0.1;                %m^2/s

%Model Start Parameters
cagesPerTreatment = 1;
noTreatments = 12;
treatmentInterval = 24;
TreatmentDepth = TreatmentDepth;

%Azamethiphos
treatmentConc = 100;    %ug/l
halfLife = 8.6;         %days
EQS72hrs = 0.04;        %ug/l
EQSPatchArea = 0.04;    %ug/l
Time = 84;              %hrs

modelTstep = 1;         %hrs
Ts = (noTreatments*treatmentInterval+Time)*3600;
T = 12.42*3600;
Phase = 0;
t=0:modelTstep*3600:Ts;

%Model Grid
gridLength = 30*1000;
gridWidth = 30*1000; 
dx = 100;
dy = 100;
X = -gridLength/2+(dx*0.5):dx:gridLength/2-(dx*0.5);
Y = -gridWidth/2+(dy*0.5):dy:gridWidth/2-(dy*0.5);
%X = 0:dx:gridLength-(dx*0.5);
%Y = 0:dy:gridWidth-(dy*0.5);
[Xgrid, Ygrid] = meshgrid(X,Y);

%Treatment
treatmentVol = CageVolume*NoCages;
TotalMass = treatmentVol*treatmentConc*1e-6;
MassPerCage = TotalMass/NoCages;   
TreatMass = MassPerCage*cagesPerTreatment;
halfLifeT = halfLife*24*3600;
kDecay = -1*log(0.5)/halfLifeT;
%Patch location
Xi=0;
Yi=0;
%Xi = gridLength/2;
%Yi = gridWidth/2;
C = zeros(length(Y),length(X),length(t));
Ci = zeros(length(Y),length(X),length(t));
%mRelease = TreatMass*noTreatments;
mRelease = TreatMass;
idx = 1:treatmentInterval:noTreatments*treatmentInterval;
idx = idx+1;
tRelease = zeros(length(t),1);
tRelease(idx) =  mRelease;


for i = 1:noTreatments
    iRelease = i;
    TreatmentIdx = find(tRelease>0);
    startIdx = TreatmentIdx(iRelease);
    Xi=0; Yi=0;
    for n=1:length(t)
        
        age = t(n)-t(startIdx)+modelTstep*3600;
        %age = t(n)-t(startIdx);
%         if age < 1  
%             age = 0;
%         end
        
        Xi = Xi+(Upara+Uamp*sin(2*pi*(t(n)/T)+Phase))*age;
        Yi = Yi+(Vpara+Vamp*sin(2*pi*(t(n)/T)+Phase))*age;
        
        Ci(:,:,n,i) = (mRelease)/(4*pi*mixingDepth*D*age)*exp(-(abs(Xgrid-Xi)/(4*D*age))-(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
        max(max(Ci(:,:,n,i)))
        
    end
end







% for i = 1:noTreatments
%     iRelease = i;
%     TreatmentIdx = find(tRelease>0);
%     startIdx = TreatmentIdx(iRelease);
%     Xi=0; Yi=0;
%     for n=1:length(t)
%         
%         %age = t(n)-t(startIdx)+modelTstep*3600;
%         age = t(n)-t(startIdx);
%         if age < 3600
%              age = 3600;
%         end
%         Xi = Xi+(Upara+Uamp*sin(2*pi*(t(n)/T)+Phase))*modelTstep*3600;
%         Yi = Yi+(Vpara+Vamp*sin(2*pi*(t(n)/T)+Phase))*modelTstep*3600;
%         
%         Ci(:,:,n,i) = (mRelease)/(4*pi*mixingDepth*D*age)*exp(-(abs(Xgrid-Xi)/(4*D*age))-(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
% max(max(Ci(:,:,n,i)))
%         
%         
%         
% %         if age < 3600
% %             age = 3600;
% %             Xi = Xi+(Upara+Uamp*sin(2*pi*(t(n)/T)+Phase))*modelTstep*3600;
% %             Yi = Yi+(Vpara+Vamp*sin(2*pi*(t(n)/T)+Phase))*modelTstep*3600;
% %             %Ci(:,:,n,i)= zeros(length(Y),length(X));
% %             Ci(:,:,n,i) = (mRelease)/(4*pi*mixingDepth*D*age)*exp(-(abs(Xgrid-Xi)/(4*D*age))-(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
% %         else
% %             Xi = Xi+(Upara+Uamp*sin(2*pi*(t(n)/T)+Phase))*modelTstep*3600;
% %             Yi = Yi+(Vpara+Vamp*sin(2*pi*(t(n)/T)+Phase))*modelTstep*3600;
% %             
% %             %Ci(:,:,n,i) = (tRelease(1)*1e+6)/(4*pi*mixingDepth*D*age)*exp(-(abs(Xgrid-Xi)/(4*D*age))-(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
% %             %Ci(:,:,n,i) = (tRelease(1)*1e+6)/(4*pi*mixingDepth*D*age)*exp(-(abs(Xgrid-Xi)/(4*D*age))-(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
% %            
% %             Ci(:,:,n,i) = (mRelease)/(4*pi*mixingDepth*D*age)*exp(-(abs(Xgrid-Xi)/(4*D*age))-(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
% %             
% %             Ci(:,:,n,i) = (mRelease)/(4*pi*mixingDepth*D*age)*exp((abs(Xgrid-Xi)/(4*D*age))+(abs(Ygrid-Yi)/(4*D*age))-kDecay*age);
% %                       
% %         end
%     end
% end



C = sum(Ci,4);
peakConc = reshape(max(C,[],[1 2]),size(C,3),1); %kg/l
%peakConc = reshape(max(Ci,[],[1 2 4]),size(Ci,3),1); %ug/l
peakConc = peakConc*1e6; %ug/l

chemMass = reshape(sum(C,[1 2]),size(C,3),1);%*1e-6; %kg


%EQS
level = 40/1000; %ng/
%Area
CLIM = 0.04;
% level = CLIM*1e-6/10;
% level = CLIM*1e-6;
for ii=1:length(t)
Area40(ii,1) = numel(find(C(:,:,ii) > level))*dx*dy/1e6;
end

[peakConc Area40 chemMass]

figure;hold on
surf(Xgrid,Ygrid,C(:,:,ii),'edgecolor','none')
view(2);caxis([0 3.7]);
contour(Xgrid,Ygrid,C(:,:,ii),[level level])



% for ii=1:length(t)
%     c = contourc(X, Y, C(:,:,ii),[level level]);
%     polygons = struct('XData', [], 'YData', []);
%     i = 1;
%     n = 0;
%     while i <= size(c,2)
%         n = n + 1;
%         l = c(2,i);
%         polygons(n).XData = c(1,i+1:i+l);
%         polygons(n).YData = c(2,i+1:i+l);      
%         i = i + l + 1;
%     end
%     a = 0;
%     for i = 1:n
%         thisArea = polyarea(polygons(i).XData, polygons(i).YData);
%         a = a + thisArea;
%     end
%     Area40(ii,1) = a/1e6;
% end

Patch = table();
Patch.Tstep = (1:length(t))';
Patch.Time = t'/(24*3600);
Patch.peakConc=peakConc;
Patch.Area40=Area40;
Patch.chemMass=chemMass;

figure
plot(peakConc)
hold on 
plot(p(1:end,1))

figure
hold on 
plot(p(1:end,1))
plot(reshape(max(Ci(:,:,:,1),[],[1 2]),size(Ci,3),1)*1e6)
plot(reshape(max(Ci(:,:,:,2),[],[1 2]),size(Ci,3),1)*1e6)
plot(reshape(max(Ci(:,:,:,3),[],[1 2]),size(Ci,3),1)*1e6)
plot(reshape(max(Ci(:,:,:,4),[],[1 2]),size(Ci,3),1)*1e6)
plot(reshape(max(Ci(:,:,:,5),[],[1 2]),size(Ci,3),1)*1e6)



figure
plot(Area40)
hold on 
plot(p(2:end,2))

figure
plot(chemMass)
hold on 
plot(p(2:end,3))


