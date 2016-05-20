#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "auxiliary.c"
#include "lmoments.c"
#include "pdfs.c"
#include "thornthwaite.c"

// Max size of raw rainfall and events matrices
#define NUMDATOSMAX 5000
#define NUMRESULTMAX 5000
#define NUMSEASONSMAX 12

// Define max() and min() functions
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

// Function prototypes
void spei(float dataSeries[], int n, int seasons, float speiSeries[]);

// Main program:
// Calculate the Standardized Precipitation Index
//int main(int argc, char **argv) {
int main(int argc, char **argv) {

	FILE *entrada,*salida;
	char pathOrigen[30],pathDestino[30],estacion[36],latitud[6];
	float lat,rainSeries[NUMDATOSMAX],tempSeries[NUMDATOSMAX],
		  etpSeries[NUMDATOSMAX],balanceSeries[NUMDATOSMAX],
		  acumSeries[NUMDATOSMAX],speiSeries[NUMDATOSMAX];
	int anio,mes,seasonality,acumulated,numRegistros,acumRegistros,
		indice,jndice;

	// Initialize variables
	anio=mes=seasonality=acumulated=numRegistros=acumRegistros=indice=jndice=0;
	lat=0.0;
	for (indice=0; indice<NUMDATOSMAX; indice++) {
		rainSeries[indice]=tempSeries[indice]=etpSeries[indice]=
			balanceSeries[indice]=acumSeries[indice]=speiSeries[indice]=0.0;
	}

	// Read in-line arguments
	if (argc!=4) {
		printf("\nUsage:\tspei [acumulated] [source file] [result file]\n");
		exit(1);
	}
	sscanf(argv[1], "%u", &acumulated);
	sscanf(argv[2], "%s", pathOrigen);
	sscanf(argv[3], "%s", pathDestino);

	// Open input file
	if((entrada=fopen(pathOrigen,"rt"))==NULL)
		{
		printf("\nError: File can't be opened");
		exit(1);
		}
	// Read heading
	fgets(estacion, 36, entrada);
	fscanf(entrada, "%f\n", &lat);
	fscanf(entrada, "%u;%u\n", &anio,&mes);
	fscanf(entrada, "%u\n", &seasonality);
	if(seasonality>NUMSEASONSMAX) {
		printf("\nError: Too many seasons. Maximum is %d", NUMSEASONSMAX);
		exit(1);
	}
	// Read data
	indice=0;
	while(!feof(entrada)) {
		if(indice==NUMDATOSMAX) {
			printf("\nError: Too many data in input file. Maximum is %d", NUMDATOSMAX);
			exit(1);
		}
		fscanf(entrada,"%f;%f\n", &rainSeries[indice],&tempSeries[indice]);
		indice++;
	}
	numRegistros=indice;
	if(tempSeries[1]==0) numRegistros-=1;
	// Close file
	fclose(entrada);
	// Print metadata (just to check)
	printf("\nseries: %s", estacion);
	printf("latitude: %.3f\n", lat);
	printf("initial date: %d/%d\n", mes, anio);
	printf("seasonality: %d\n", seasonality);
	printf("%d registers\n", numRegistros);
	printf("calculating SPEI at %d month", acumulated);
	if (acumulated>1) printf("s");
	printf("\n");

	// Compute the climatic balance: precipitation minus potential evapotranspiration
	if (tempSeries[1]!=0 && tempSeries[2]!=0) {
		thornthwaite(tempSeries, numRegistros, lat, etpSeries);
		for (indice=0; indice<numRegistros; indice++) {
			balanceSeries[indice] = rainSeries[indice]-etpSeries[indice];
		}
	}
	else {
		for (indice=0; indice<numRegistros; indice++) {
			balanceSeries[indice] = rainSeries[indice];
		}
	}

	// Compute the cumulative series
	anio += (acumulated-1)/12;
	mes += acumulated-1;
	while (mes>12) mes-=12;
	acumRegistros = numRegistros-acumulated+1;
	for (indice=acumulated-1; indice<numRegistros; indice++) {
		for (jndice=0; jndice<acumulated; jndice++) {
			acumSeries[indice-acumulated+1] += balanceSeries[indice-jndice];
		}
	}

	// Compute the SPEI series
	spei(acumSeries, acumRegistros, seasonality, speiSeries);

	// Write results to file
	if((salida=fopen(pathDestino,"wt"))==NULL) {
		printf("\nError: Output file could not be opened");
		exit(1);
	}
	fprintf(salida,"%s%f\n%u;%u\n%u", estacion,lat,anio,mes,seasonality);
	//for (jndice=1; jndice<=seasonality; jndice++) {
	//	fprintf(salida,"\n%f;%f;%f", logLogisticParams[jndice][0],
	//			logLogisticParams[jndice][1], logLogisticParams[jndice][2]);
	//}
	for (indice=0; indice<acumRegistros; indice++) {
		fprintf(salida,"\n%f", speiSeries[indice]);
	}
	fclose(salida);

	// Quit
	return(0);
}

// spei()
// Calculates the Standardized Precipitation-Evapotransporation Index
// from a series of climatic balance (precipitation minus etp). The
// SPEI is the standardized value of the climatic balance (P-ETP),
// computed following a Log Logistic probability distribution.
void spei(float dataSeries[], int n, int seasons, float speiSeries[]) {

	int i, j, k, nSeason;
	float seasonSeries[NUMDATOSMAX], beta[3], logLogisticParams[NUMSEASONSMAX][3];

	// Loop through all seasons defined by seasons
	for (j=1; j<=seasons; j++) {
		// Extract and sort the seasonal series
		k = 0;
		for (i=j-1; i<n; i+=seasons) {
			seasonSeries[k] = dataSeries[i];
			k++;
		}
		nSeason = k;
		upward(seasonSeries, nSeason);
		// Compute probability weighted moments
		//pwm(seasonSeries, nSeason, beta, -0.35, 0, 0);
		pwm(seasonSeries, nSeason, beta, 0, 0, 0);
		// Fit a Log Logistic probability function
		logLogisticFit(beta, logLogisticParams[j]);
		//printf("\nSeason %u", jndice);
		//printf("\nLogLogistic beta param.: %.4f", logLogisticParams[jndice][0]);
		//printf("\nLogLogistic alpha param.: %.4f", logLogisticParams[jndice][1]);
		//printf("\nLogLogistic gamma param.: %.4f\n", logLogisticParams[jndice][2]);
		// Calculate the standardized values
		for (i=j-1; i<n; i+=seasons) {
			speiSeries[i] = logLogisticCDF(dataSeries[i], logLogisticParams[j]);
			speiSeries[i] = -standardGaussianInvCDF(speiSeries[i]);
		}
	}
}
