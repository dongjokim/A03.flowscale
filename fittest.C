#include "vnPtIntegrate.C"

void fittest(){
	//for(int ic=0;ic<NC;ic++) {
	for(int ic=0;ic<1;ic++) {
		FitPtMinuit(ic);
	}
}
