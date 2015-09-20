#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Function.h"

#define NAME_IMG_IN  "photograph"
#define NAME_IMG_OUT1 "photograph_noise"
#define NAME_IMG_OUT2 "photograph_denoise"


int main(int argc,char** argv){

  	int i,j;
  	int length,width;
  	int nbLevels;
	
	float** image;
	float** haar;
	float** temp;
	float** temp1;
	float** haar_inverse;
	float threshold;
	float sum1,sum2;
	float var;   
	
	printf("Input the variance of gaussian noise: ");
  	scanf("%f",&var);
	
	printf("Input the number of level: ");
  	scanf("%d",&nbLevels);
		
	printf("Input the threshold: ");
  	scanf("%f",&threshold);

	image=LoadImagePgm(NAME_IMG_IN, &length, &width);
	haar=fmatrix_allocate_2d(length,width);
	temp=fmatrix_allocate_2d(length,width);
	temp1=fmatrix_allocate_2d(length,width);
	haar_inverse=fmatrix_allocate_2d(length,width);

	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			haar[i][j]=0.0;
			temp[i][j]=0.0;
			temp1[i][j]=0.0;
			haar_inverse[i][j]=0.0;
		}
	}
	/* add_gaussian_noise */
	copy(temp,image,length,width);
	add_gaussian_noise(temp,length,width,var);

 	/* denoise image based on Haar coefficients */
	haar2D_complete(temp,haar,nbLevels,length,width);
	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			if((i>length/pow(2,nbLevels)||j>width/pow(2,nbLevels))){
				if(fabs(haar[i][j])<threshold){
					haar[i][j]=0.0;
				}
				else{
					haar[i][j]=255*(haar[i][j]/fabs(haar[i][j]))*((fabs(haar[i][j])-threshold)/(255-threshold));
				}
			}
			else{
				haar[i][j]=haar[i][j];
						
			}		
		}
	}

	ihaar2D_complete(haar,haar_inverse,nbLevels,length,width);	
	Recal(haar_inverse,length,width);

	/* display ISNR */
	sum1=0.0;
	sum2=0.0;
	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			sum1=sum1+SQUARE(image[i][j]-temp[i][j]);
			sum2=sum2+SQUARE(image[i][j]-haar_inverse[i][j]);
		}
	}
	
	printf("\nISNR=%f\n",10*log10(sum1/sum2));
 
 	/* save images */
 	SaveImagePgm(NAME_IMG_OUT1,temp,length,width);
	system("display photograph_noise.pgm&");
	SaveImagePgm(NAME_IMG_OUT2,haar_inverse,length,width);
	system("display photograph_denoise.pgm&");

	/*free memory*/
	free_fmatrix_2d(image);  
	free_fmatrix_2d(haar);
	free_fmatrix_2d(temp);
	free_fmatrix_2d(temp1);
	free_fmatrix_2d(haar_inverse);
  	
  	printf("\n Ending ... \n\n\n");
  	return 0; 	 
}
