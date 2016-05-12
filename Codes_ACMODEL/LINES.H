#ifndef _LINES
#define _LINES

void SuctionLine(HP,double,HP*,MASS*);
void SuctionLine_Lumped(HP,double,HP*,MASS*);
void DistributionTubes(HP,double,HP*,MASS*);
void HotGasLine(HP,double,HP*,MASS*,void*);
void HotGasLine_Lumped(HP,double,HP*,MASS*,void*);
void LiquidLine(HP,double,HP*,MASS*);
void LiquidLine_Lumped(HP,double,HP*,MASS*);
void Distribution_CapTubes(HP,double,HP*,MASS*,double,double,int);
#endif
