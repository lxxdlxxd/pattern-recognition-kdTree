//*********************************************************Oct 6th 2017, version 1.0, author Xin Li***************************************************************
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<algorithm>
#include<limits>
using namespace std;

#define NUMMDRF 252 // 81 //252 
#define NUMPMT 20 // 16 //20 // size of the data vector, determined by number of photosensors

class vectorK{ // the data structure to store both the index/value of a specific reference data vector
public:
    int idx;
    float value[NUMPMT];
};

class node{ // the node data structure, idx indicates the axis to split, leaf indicates wether it is a internal node "-1" or leaf node, it also stores the index of the reference data vector in the MDRFs if it is a leaf node.
public:
    int idx, leaf;
    float boundary;
    node *leftChild;
    node *rightChild;
};

void constructTree(vectorK *vectors, node *root, int size){ // the function to construct the kd search tree, "vectors" contains the reference data to be stored in the tree, "root" is the pointer to the tree, size is the size of "vectors"
    if(size>1){ //If the size of "vectors" is larger than 1, then it is an internal node. The internal node records the axis (x or y or z, etc) and boundary position (x=0.3 or y=0.25 or z=0.7, etc) to split the data.
        int i, j, n=size/2, record=0, idxChildLeft=0, idxChildRight=0, iter, countSameLeft, countSameLeftRecord;
        float boundary, temp, maxi=numeric_limits<float>::min(), boundaryRecord;
        vector<float> vectorTemp(size);
        vectorK *vectorChildLeft, *vectorChildRight; //vectorChildLeft stores the data on the left side of the boundary, vectorChildRight stores the data on the right side of the boundary.
        vectorChildLeft=(vectorK *)malloc(sizeof(vectorK)*n);
        vectorChildRight=(vectorK *)malloc(sizeof(vectorK)*(size-n));
        for(i=0;i<NUMPMT;i++){
            for(j=0;j<size;j++){
                vectorTemp[j]=vectors[j].value[i]; //"vectorTemp" stores the ith element of PMT values for all the events in "vectors"
            }
            std::sort(vectorTemp.begin(), vectorTemp.end()); //sort the elements in "vectorTemp"
            boundary = (vectorTemp[n]+vectorTemp[n-1])/2; // find the boundary to split the data, n indicates the median of the "vectorTemp"
            countSameLeft=0;
            if(vectorTemp[n-1]==vectorTemp[n]){
                iter=n-1;
                while(iter>=0&&vectorTemp[iter]==vectorTemp[n-1]){ //continues to count the same elements until the element is different.
                    iter--;
                    countSameLeft++; //countSameLeft is used to balance the branches on the left and right.
                }
            }
            temp=0;
            for(j=0;j<size;j++){
                temp+=pow((vectorTemp[j]-boundary), 2); // a metric used to find the axis with the largest split.
            }
            if(temp>maxi){
                record=i; // record the axis index if it is the current largest.
                maxi=temp;
                boundaryRecord=boundary; // record the boundary
                countSameLeftRecord=countSameLeft; // record the countSameLeft
            }
        }
        root->idx=record; // store the axis index into the kd tree.
        root->boundary=boundaryRecord; // store the boundary into the kd tree.
        for(i=0;i<size;i++){
            if(vectors[i].value[record]<boundaryRecord){
                vectorChildLeft[idxChildLeft].idx=vectors[i].idx; // pass the location of reference vector in the MDRF into "vectorChildLeft", which is variable "vectors" for the left child or left branch.
                for(j=0;j<NUMPMT;j++){
                    vectorChildLeft[idxChildLeft].value[j]=vectors[i].value[j]; // pass the values of the reference vector into its left child or left branch.
                }
                idxChildLeft++;
            }else if(vectors[i].value[record]>boundaryRecord){
                vectorChildRight[idxChildRight].idx=vectors[i].idx;
                for(j=0;j<NUMPMT;j++){
                    vectorChildRight[idxChildRight].value[j]=vectors[i].value[j];
                }
                idxChildRight++;
            }else{ // if the data vector is on the boundary, then assign it to either left or right boundary to balance the left and right branch.
                if(countSameLeftRecord){ 
                    vectorChildLeft[idxChildLeft].idx=vectors[i].idx;
                    for(j=0;j<NUMPMT;j++){
                        vectorChildLeft[idxChildLeft].value[j]=vectors[i].value[j];
                    }
                    idxChildLeft++;
                    countSameLeftRecord--;
                }else{
                    vectorChildRight[idxChildRight].idx=vectors[i].idx;
                    for(j=0;j<NUMPMT;j++){
                        vectorChildRight[idxChildRight].value[j]=vectors[i].value[j];
                    }
                    idxChildRight++;
                }
            }
        }
        vectorTemp.clear(); // destroy the "vectorTemp" variable
        vector<float>().swap(vectorTemp); // further clean the "vectorTemp"
        delete [] vectors; // free(vectors);
        root->leftChild=(node *)malloc(sizeof(node));
        root->leftChild->leaf=-1; // if it is internal node, mark it as -1
        constructTree(vectorChildLeft, root->leftChild, n);
        root->rightChild=(node *)malloc(sizeof(node));
        root->rightChild->leaf=-1;
        constructTree(vectorChildRight, root->rightChild, size-n); // contruct kd tree iteratively
    }else{ //If the size is 1 or 0, the node is a leaf node, the leaf node record the index "idx" or location of the data vector in the MDRF.
        root->leaf=vectors->idx; // if it is leaf node, record the location on the MDRF.
        delete [] vectors;
    }
}

void saveTree(node *root, FILE *fid){ // save the kd tree to disk
    fwrite(&(root->idx), sizeof(int), 1, fid);
    fwrite(&(root->boundary), sizeof(float), 1, fid);
    fwrite(&(root->leaf), sizeof(int), 1, fid);
    if(root->leaf==-1){
        saveTree(root->leftChild, fid); // save the kd tree iteratively
        saveTree(root->rightChild, fid);
    }
}

void MDRF_filteration(float *MDRF_filtered, float *MDRF){ // filter the MDRFs with average filter, but for this work it will do nothing, MDRF_filtered = MDRF just ignore it.
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

void normalize(float *MDRF_filtered_try){// normalize the MDRFs
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

int main(){
    int i, j, size = NUMMDRF*NUMMDRF;
    float *MDRF, *MDRFRead;
    FILE *fid;
    vectorK *vectors;
    node *root;

    MDRF=(float *)malloc(sizeof(float)*NUMMDRF*NUMMDRF*NUMPMT);
    MDRFRead=(float *)malloc(sizeof(float)*NUMMDRF*NUMMDRF*NUMPMT);
    vectors=(vectorK *)malloc(sizeof(vectorK)*NUMMDRF*NUMMDRF);
    root = (node *)malloc(sizeof(node));
    root->leaf=-1;

    fid = fopen("E:\\research\\data\\20180314_cnn\\MDRFSimulated_252X252_NumPhoton200000_NumSiPMs20_CsI_Barrier1_QE0.2_lambertianRef0.9_gelIndex1.500_WholeLargeRodSpace1_.bin", "rb"); //E:\\research\\data\\20180115\\MDRFAnnealed5_50000_50000_5000000.bin
    if(fread(MDRFRead, sizeof(float), NUMMDRF*NUMMDRF*NUMPMT, fid)){
        cout<<"Successfully read MDRF"<<endl;
    }else{
        cout<<"Error reading MDRF"<<endl;
        exit(0);
    }
    fclose(fid);

    MDRF_filteration(MDRF, MDRFRead);
	normalize(MDRF);

    for(i=0;i<NUMMDRF*NUMMDRF;i++){
        vectors[i].idx=i;
        for(j=0;j<NUMPMT;j++){           
            vectors[i].value[j]=MDRF[i*NUMPMT+j];
        }
    }

    constructTree(vectors, root, size);
    fid=fopen("E:\\research\\data\\20180508_7lineTrueLable\\treeSim.bin", "wb");
    saveTree(root, fid);
    fclose(fid);
    return 0;
}