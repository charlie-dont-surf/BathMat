# BathMat
BathMat is a tool used to calculate maximum chemical treatment quantities using basic assumptions for short term EQS. The software uses identical methods as that implemented in BathAuto. The results have been compared and BathMat produces identical results. 


To use BathMat

1. Open BathMat.m and input site cage and flow parameters. For the short-term model only Umean is required. Strictly speaking SiteName and SiteCenter are not required to run the model.

2. Apply short term EQS values. Time must be in hours. EQSconc and TreatmentConc must be in ng/l.
 
3. Use function BathMat_ShortTermModel to calculate maximum chemical quantities


Additional code is provided below for a Long-term Model (Gillibrand and Turrell 1999). However, this is not functional at this time. 
