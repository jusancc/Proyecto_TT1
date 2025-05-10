#include "..\include\global.hpp"


/*extern*/ Matrix eopdata;

void eop19620101(int c){
	eopdata = zeros(13, c);

	FILE *fid = fopen("../data/eop19620101.txt", "r");
	if(fid == NULL){
		printf("Fail open o file\neop19620101");
		exit(EXIT_FAILURE);
	}

	for (int j=1; j<=c; j++){
		fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(eopdata(1,j)),&(eopdata(2,j)),&(eopdata(3,j)),
		&(eopdata(4,j)),&(eopdata(5,j)),&(eopdata(6,j)),&(eopdata(7,j)),&(eopdata(8,j)),&(eopdata(9,j)),&(eopdata(10,j)),
		&(eopdata(11,j)),&(eopdata(12,j)),&(eopdata(13,j)));
	}
	
}

Matrix Cnm, Snm;
void GGM03S(int c){
	Cnm = zeros(181, 181);
	Snm = zeros(181, 181);
	double aux;

	FILE *fid = fopen("../data/GGm03S.txt", "r");
	for (int i=0;i<180;i++){
		for (int j=0;j<i;j++){
			fscanf(fid, "%lf %lf %lf %lf %lf %lf", &aux, &aux, Cnm(i+1,j+1), Snm(i+1,j+1), &aux, &aux);
		}
	}
}

Matrix PC;
void DE430Coeff(int f, int c){
    PC = zeros(f,c);

    FILE *fid = fopen("../data/DE430Coeff.txt","r");
    if(fid == NULL){
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    for(int i=1; i<=f; i++){
        for(int j=1; j<=c; j++){
            fscanf(fid,"%lf",&(PC(i,j)));
        }
    }

    fclose(fid);
}