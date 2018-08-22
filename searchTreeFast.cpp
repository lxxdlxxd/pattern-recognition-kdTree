//*********************************************************Oct 6th 2017, version 1.0, author Xin Li***************************************************************
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<algorithm>
#include<limits>
#include<ctime>
using namespace std;

#define NUMMDRF 252 // 81 //81 //252
#define NUMPMT 20 // 16 //16 //20

#define INTERVALMDRF 0.2 // 0.25 // 0.2 //0.25

class node{
public:
    int idx, leaf;
    float boundary;
    node *leftChild;
    node *rightChild;
};

int record;
float *MDRF, budgetTotal, halfMDRF=(NUMMDRF-1)*INTERVALMDRF/2;

void loadTree(node *root, FILE *fid){ // read the tree from disk, and load it into variable "root"
    fread(&(root->idx), sizeof(int), 1, fid);
    fread(&(root->boundary), sizeof(float), 1, fid);
    fread(&(root->leaf), sizeof(int), 1, fid);
    if(root->leaf==-1){
        root->leftChild=(node *)malloc(sizeof(node));
        loadTree(root->leftChild, fid);
        root->rightChild=(node *)malloc(sizeof(node));
        loadTree(root->rightChild, fid);
    }
}

void MDRF_filteration(float *MDRF_filtered, float *MDRF){ // I do nothing, just ignore me.
	int idxPMT, i, j, flag;

	for (idxPMT = 0; idxPMT < NUMPMT; idxPMT++){
		for (i = 0; i < NUMMDRF; i++){
			for (j = 0; j < NUMMDRF; j++){
				MDRF_filtered[(i*NUMMDRF + j)*NUMPMT + idxPMT] = MDRF[(i*NUMMDRF + j)*NUMPMT + idxPMT];
			}
		}
	}
	for (idxPMT = 0; idxPMT < NUMPMT; idxPMT++){
		for (i = 1; i < NUMMDRF - 1; i++){
			for (j = 1; j < NUMMDRF - 1; j++){
				flag = MDRF[(i*NUMMDRF + j)*NUMPMT + idxPMT] * MDRF[((i - 1)*NUMMDRF + j)*NUMPMT + idxPMT] * MDRF[(i*NUMMDRF + j - 1)*NUMPMT + idxPMT] * MDRF[((i + 1)*NUMMDRF + j)*NUMPMT + idxPMT] * MDRF[(i*NUMMDRF + j + 1)*NUMPMT + idxPMT] * MDRF[((i - 1)*NUMMDRF + j - 1)*NUMPMT + idxPMT] * MDRF[((i - 1)*NUMMDRF + j + 1)*NUMPMT + idxPMT] * MDRF[((i + 1)*NUMMDRF + j - 1)*NUMPMT + idxPMT] * MDRF[((i + 1)*NUMMDRF + j + 1)*NUMPMT + idxPMT] != 0;
				if (flag == 0){
					MDRF_filtered[(i*NUMMDRF + j)*NUMPMT + idxPMT] = MDRF[(i*NUMMDRF + j)*NUMPMT + idxPMT];
				}
				else{
					MDRF_filtered[(i*NUMMDRF + j)*NUMPMT + idxPMT] = MDRF[(i*NUMMDRF + j)*NUMPMT + idxPMT]; //(MDRF[(i*NUMMDRF + j)*NUMPMT + idxPMT] * 6 + MDRF[((i - 1)*NUMMDRF + j)*NUMPMT + idxPMT] * 2 + MDRF[(i*NUMMDRF + j - 1)*NUMPMT + idxPMT] * 2 + MDRF[((i + 1)*NUMMDRF + j)*NUMPMT + idxPMT] * 2 + MDRF[(i*NUMMDRF + j + 1)*NUMPMT + idxPMT] * 2 + MDRF[((i - 1)*NUMMDRF + j - 1)*NUMPMT + idxPMT] + MDRF[((i - 1)*NUMMDRF + j + 1)*NUMPMT + idxPMT] + MDRF[((i + 1)*NUMMDRF + j - 1)*NUMPMT + idxPMT] + MDRF[((i + 1)*NUMMDRF + j + 1)*NUMPMT + idxPMT]) / 18.0; //MDRF[((i + 1)*NUMY + j + 1)*NUMPMT + idxPMT] + MDRF[((i - 1)*NUMY + j + 1)*NUMPMT + idxPMT] + MDRF[((i - 1)*NUMY + j - 1)*NUMPMT + idxPMT] + MDRF[((i-1)*NUMY + j - 1)*NUMPMT + idxPMT]
				}
			}
		}
	}
}

void normalize(float *MDRF_filtered_try){
	int i, j, idxPMT;
	float temp;

	for (i = 0; i < NUMMDRF; i++){
		for (j = 0; j < NUMMDRF; j++){
			temp = 0;
			for (idxPMT = 0; idxPMT < NUMPMT; idxPMT++){
				temp = temp + MDRF_filtered_try[(i*NUMMDRF + j)*NUMPMT + idxPMT];
			}

			for (idxPMT = 0; idxPMT < NUMPMT; idxPMT++){
				MDRF_filtered_try[(i*NUMMDRF + j)*NUMPMT + idxPMT] = MDRF_filtered_try[(i*NUMMDRF + j)*NUMPMT + idxPMT] / temp;
			}
		}
	}
}

void searchTree(node *root, float *event, float *accumulateCostParent, float costTotalSquareParent){
    int idx=root->idx;
    if(root->leaf==-1){ // judge whether it is a internal node
        float boundary=root->boundary, accumulateCost[NUMPMT], costTotalSquare;
		memcpy(accumulateCost, accumulateCostParent, sizeof(float)*NUMPMT);
		costTotalSquare = costTotalSquareParent;
        if(event[idx]<boundary){
            float test=event[idx]; // a variable for testing, just ignore
            searchTree(root->leftChild, event, accumulateCost, costTotalSquare); // search the left branch
			if (accumulateCost[idx] < boundary - event[idx]) {
				costTotalSquare -= pow(accumulateCost[idx], 2);
				accumulateCost[idx] = boundary - event[idx];
				costTotalSquare += pow(accumulateCost[idx], 2);
			}
            if(budgetTotal>costTotalSquare){ // if the accumulated error "budgetTotal" is larger than the euclidean distance to the boundary, also search the other branch
                searchTree(root->rightChild, event, accumulateCost, costTotalSquare);
            }
        }else{
            float test=event[idx]; // a variable for testing, just ignore
            searchTree(root->rightChild, event, accumulateCost, costTotalSquare); // search the right branch
			if (accumulateCost[idx] < event[idx] - boundary) {
				costTotalSquare -= pow(accumulateCost[idx], 2);
				accumulateCost[idx] = event[idx] - boundary;
				costTotalSquare += pow(accumulateCost[idx], 2);
			}
            if(budgetTotal>costTotalSquare){ // if the accumulated error "budgetTotal" is larger than the euclidean distance to the boudary, also search the other branch
                searchTree(root->leftChild, event, accumulateCost, costTotalSquare);
            }
        }
    }else{
        int i, idxMDRF=root->leaf;
        float temp=0;

        for(i=0;i<NUMPMT;i++){
            temp+=pow(MDRF[idxMDRF*NUMPMT+i]-event[i], 2); // calculate the accumulated error
        }

        if(temp<budgetTotal){ // update the accumulated error with the minimum error found so far
            budgetTotal=temp;
            record=idxMDRF; // record the index of the reference vector who has minimum error so far
        }
    }
}

int main(){
    int flag, size, i, j;
    float x, y, *MDRFRead, accumulateCost[NUMPMT], costTotalSquare, *PMTSignals, *event, total;
    node *root;
    FILE *fid;

    MDRF=(float *)malloc(sizeof(float)*NUMMDRF*NUMMDRF*NUMPMT);
    MDRFRead=(float *)malloc(sizeof(float)*NUMMDRF*NUMMDRF*NUMPMT);
    event=(float *)malloc(sizeof(float)*NUMPMT);

    root=(node *)malloc(sizeof(node));

    fid=fopen("E:\\research\\data\\20180508_7lineTrueLable\\treeSim.bin", "rb");
    loadTree(root, fid); // load the tree
    fclose(fid);

    fid = fopen("E:\\research\\data\\20180314_cnn\\MDRFSimulated_252X252_NumPhoton200000_NumSiPMs20_CsI_Barrier1_QE0.2_lambertianRef0.9_gelIndex1.500_WholeLargeRodSpace1_.bin", "rb"); //E:\\research\\data\\20180115\\MDRFAnnealed5_50000_50000_5000000.bin
	flag = fread(MDRFRead, sizeof(float), NUMMDRF*NUMMDRF*NUMPMT, fid);
    if (!flag){
		printf("Cannot read MDRF.\n");
		exit(0);
	}
	else{
		printf("MDRF loaded successfully.\n");
	}
    fclose(fid);
    
    MDRF_filteration(MDRF, MDRFRead);
    normalize(MDRF);
    
    fid = fopen("E:\\research\\data\\20180205\\7lineSimulated.bin", "rb"); // the test data, detector responses by illuminating the detector with 5 gamma-ray slit beams. //E:\\research\\data\\20180109\\calibrationV_73.bin
    fseek(fid, 0, SEEK_END);
	size = ftell(fid)/sizeof(float)/NUMPMT;
	PMTSignals = (float *)malloc(sizeof(float)*size*NUMPMT);
	fseek(fid, 0, SEEK_SET);
    flag = fread(PMTSignals, sizeof(float), size*NUMPMT, fid);
    if (!flag){
		printf("Cannot read PMTSignals.\n");
		exit(0);
	}
	else{
		printf("PMTSignals loaded successfully.\n");
	}
    fclose(fid);
    
    fid = fopen("E:\\research\\data\\20180508_7lineTrueLable\\resultKDTreeFast.bin", "wb"); //E:\\research\\data\\20180205\\resultKDTreeFast.bin

    std::clock_t start, end;
    start = std::clock();
    
    for(i=0;i<size;i++){
        if(i%100000==0){
            cout<<i<<" is finished."<<endl;
        }
		memset(accumulateCost, 0, sizeof(float)*NUMPMT);
		costTotalSquare = 0;
        budgetTotal=numeric_limits<float>::max(); // initialize the accumulated error with the maximum float num
        total=0;
        for(j=0;j<NUMPMT;j++){
            event[j]=PMTSignals[i*NUMPMT+j];
            total+=event[j];
        }
        for(j=0;j<NUMPMT;j++){
            event[j]/=total;
        }
        searchTree(root, event, accumulateCost, costTotalSquare); // search for the query vector "event"
        x=(record/NUMMDRF)*INTERVALMDRF-halfMDRF; // convert the index in the MDRF into physical location in the detector 
        y=(record%NUMMDRF)*INTERVALMDRF-halfMDRF; // convert the index in the MDRF into physical location in the detector
		/*budgetTotal = 0;
		for (j = 0; j < NUMPMT; j++) {
			budgetTotal += pow(event[j] * total - MDRF[record*NUMPMT+j] * total, 2);
		}
		budgetTotal = 1.0 / budgetTotal;*/
		budgetTotal = 1.0 / budgetTotal / pow(total, 2);
        fwrite(&budgetTotal, sizeof(float), 1, fid); // record the error
        fwrite(&x, sizeof(float), 1, fid);
        fwrite(&y, sizeof(float), 1, fid);
    }

    end = std::clock();
	printf("%ld\n", end-start);
	getchar();
    return 0;
}