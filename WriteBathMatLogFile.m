function [log] = WriteBathMatLogFile(FolderPath,Logfilename,Input,BathMat1,AzaEQS,CypEQS,DelEQS)
%% Write log file
log = {};
log{1,1} = 'BathMat - Log file';

log{3,1} = 'Site Details' ;
log{4,1} = ['SiteName: ' Input.SiteName];
log{5,1} = ['Site Location: ' num2str(Input.SiteCenter)];
log{6,1} = ['Number of Cages: ' num2str(Input.NoCages)];
log{7,1} = ['Cage Circ (m): ' num2str(Input.CageCirc)];
log{8,1} = ['Site Depth (m): ' num2str(Input.SiteDepth)];
log{9,1} = ['Distance to shore (m): ' num2str(Input.Dist2Shore)];

log{11,1} = ['TreatmentDepth (m): ' num2str(Input.TreatmentDepth)];
log{12,1} = ['cageShape = ' Input.cageShape];
log{13,1} = ['wellboatCapacity (m^3): ' num2str(Input.wellboatCapacity)];
log{14,1} = ['Mean flow speed (m/s): ' num2str(Input.Umean)];

log{16,1} = ['Azamethiphos'];
log{17,1} = ['Time (hours): ' num2str(AzaEQS.Time)];
log{18,1} = ['EQSconc (ng/l): ' num2str(AzaEQS.EQSconc)];
log{19,1} = ['Treatment Conc (ng/l): ' num2str(AzaEQS.TreatmentConc)];

log{21,1} = ['Cypermethrin'];
log{22,1} = ['Time (hours): ' num2str(CypEQS.Time)];
log{23,1} = ['EQSconc (ng/l): ' num2str(CypEQS.EQSconc)];
log{24,1} = ['Treatment Conc (ng/l): ' num2str(CypEQS.TreatmentConc)];

log{26,1} = ['Deltamethrin'];
log{27,1} = ['Time (hours): ' num2str(DelEQS.Time)];
log{28,1} = ['EQSconc (ng/l): ' num2str(DelEQS.EQSconc)];
log{29,1} = ['Treatment Conc (ng/l): ' num2str(DelEQS.TreatmentConc)];

log{31,1} = ['Results...'];
StartIdx = 32; 
for n = 1:height(BathMat1)
    log{StartIdx+n-1,1} = ['Treatment type: ' BathMat1.Chemical{n,1}];
    StartIdx = StartIdx+1; 
    if contains(BathMat1.Chemical{n,1},{'wellboat','well boat','well-boat','well_boat'},'IgnoreCase',true)
        log{StartIdx+n-1,1} = ['Cage Volume (m^3): ' num2str(Input.CageVolumeWellBoat)];
    else
        log{StartIdx+n-1,1} = ['Cage Volume (m^3): ' num2str(Input.CageVolume)];  
    end
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Plume Length (m): ' num2str(BathMat1.plumeLength(n,1))];
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Plume Width (m): ' num2str(BathMat1.plumeWidth(n,1))];
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Mixing Zone Area (m^2): ' num2str(BathMat1.mixingZoneArea(n,1))];
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Treatment Volume (m^3): ' num2str(BathMat1.treatmentVol(n,1))];
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Number of Treatable Cages: ' num2str(BathMat1.noCagesTreated(n,1))];
    StartIdx = StartIdx+1;    
    log{StartIdx+n-1,1} = ['Mean Concentration (ng/l): ' num2str(BathMat1.meanConc(n,1))];
    StartIdx = StartIdx+1;   
    log{StartIdx+n-1,1} = ['Consent Mass (g): ' num2str(BathMat1.consentMass(n,1))];
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Peak Concentration (ng/l): ' num2str(BathMat1.peakConc(n,1))];
    StartIdx = StartIdx+1;
    log{StartIdx+n-1,1} = ['Area Exceeding EQS (m^2): ' num2str(BathMat1.areaExceedsEQS(n,1))];
    StartIdx = StartIdx+2;    
    
end

%write log file
cd(FolderPath);
fileID = fopen(Logfilename,'w');
fprintf(fileID,'%s\n',log{:});
fclose(fileID);
disp(['Log File Created: ' Logfilename ' ...........'])

end

