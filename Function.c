#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Function.h"

float* fmatrix_allocate_1d(int hsize)
 {
  float* matrix;

  matrix=(float*)malloc(sizeof(float)*hsize); 
  if (matrix==NULL) printf("problem of allocating memory");

  return matrix; 
 }


float** fmatrix_allocate_2d(int vsize,int hsize)
 {
  int i;
  float** matrix;
  float *imptr;

  matrix=(float**)malloc(sizeof(float*)*vsize);
  if (matrix==NULL) printf("problem of allocating memory");

  imptr=(float*)malloc(sizeof(float)*hsize*vsize);
  if (imptr==NULL) printf("problem of allocating memory");
 
  for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }


void free_fmatrix_1d(float* pmat)
{ 
  free(pmat); 
}


void free_fmatrix_2d(float** pmat)
 { 
  free(pmat[0]);
  free(pmat);
 }


float** LoadImagePgm(char* name,int *length,int *width)
{
  int i,j,k;
  unsigned char var;
  char buff[NBCHAR];
  float** mat;

  char stringTmp1[NBCHAR],stringTmp2[NBCHAR],stringTmp3[NBCHAR];
 
  int ta1,ta2,ta3;
  FILE *fic;

  /*-----load image name-----*/
  strcpy(buff,name);
  strcat(buff,".pgm");
  printf("---> Open %s",buff);

  /*----open file----*/
  fic=fopen(buff,"r");
  if (fic==NULL)
    { printf("\n- There is en error happened when open %s  -\n",buff);
      exit(-1); }

  /*--get image info--*/
  fgets(stringTmp1,100,fic);
  fgets(stringTmp2,100,fic);
  fscanf(fic,"%d %d",&ta1,&ta2);
  fscanf(fic,"%d\n",&ta3);

  /*--display info--*/
  printf("\n\n--Info--");
  printf("\n----------");
  printf("\n%s%s%d %d \n%d\n",stringTmp1,stringTmp2,ta1,ta2,ta3);

  *length=ta1;
  *width=ta2;
  mat=fmatrix_allocate_2d(*length,*width);
   
  /*--load matrix data--*/
   for(i=0;i<*length;i++)
      for(j=0;j<*width;j++)  
        { fread(&var,1,1,fic);
          mat[i][j]=var; }

  fclose(fic);

  return(mat);
}


void SaveImagePgm(char* name,float** mat,int length,int width)
{
  int i,j,k;
  char buff[NBCHAR];
  FILE* fic;
  time_t tm;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".pgm");

  /*--open file--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Problem when save %s",buff); 
          exit(-1); }
  printf("\n Save %s in format pgm\n",name);

  /*--save info--*/
  fprintf(fic,"P5");
  if (ctime(&tm)==NULL) fprintf(fic,"\n#\n");
  else fprintf(fic,"\n# IMG Module, %s",ctime(&tm));
  fprintf(fic,"%d %d",width,length);
  fprintf(fic,"\n255\n");

  /*--save data--*/
  for(i=0;i<length;i++)
      for(j=0;j<width;j++) 
        fprintf(fic,"%c",(char)mat[i][j]);
   
   fclose(fic); 
 } 

/*--------------*/
/* FOURIER -----*/
/*--------------*/
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}


void FFTDD(float** mtxR,float** mtxI,int lgth, int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocating memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1)); 

 /*fill nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*fill data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,FFT);

 /*fill data*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2];  }

 /*matrix conversion*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

    mtxR[posy][posx]=ImgFreqR[i];  
    mtxI[posy][posx]=ImgFreqI[i];}

 /*free memory*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}


void IFFTDD(float** mtxR,float**  mtxI,int lgth,int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocating memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

 /*fill nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*fill data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,IFFT);

 /*fill data*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2]; }

 /*matrix conversion*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

   mtxR[posy][posx]=ImgFreqR[i]/(wdth*lgth);  
   mtxI[posy][posx]=ImgFreqI[i]/(wdth*lgth); }

 /*free memory*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}


void Mod(float** matM,float** matR,float** matI,int lgth,int wdth)
{
 int i,j;

 /*Calcul module*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 matM[i][j]=sqrt((matR[i][j]*matR[i][j])+(matI[i][j]*matI[i][j]));
}


void Mult(float** mat,float coef,int lgth,int wdth)
{
 int i,j;

 /*multiplication*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    { mat[i][j]*=coef;
      if (mat[i][j]>GREY_LEVEL) mat[i][j]=GREY_LEVEL; }
}


void Recal(float** mat,int lgth,int wdth)
{
 int i,j;
 float max,min;

 /*Initialization*/
 min=mat[0][0];

 /*look for min*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if (mat[i][j]<min) min=mat[i][j];

 /*plus min*/
   for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][j]-=min;

   max=mat[0][0];
 /*look for max*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    if (mat[i][j]>max) max=mat[i][j];

 /*matrix caliberation*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   mat[i][j]*=(GREY_LEVEL/max);      
}


void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
     matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}


void SquareMatrix(float** matRout,float** matIout,float** matRin,float** matIin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=SQUARE(matRin[i][j])-SQUARE(matIin[i][j]);
     matIout[i][j]=2*matRin[i][j]*matIin[i][j]; }
}

/*****************************************************/
/* gaussian noise                                    */
/*****************************************************/
float gaussian_noise(float var,float mean)
{
 float noise,theta;

 /*generation of noise*/
 noise=sqrt(-2*var*log(1.0-((float)rand()/RAND_MAX)));
 theta=(float)rand()*1.9175345E-4-PI;
 noise=noise*cos(theta);
 noise+=mean;
 if (noise>GREY_LEVEL) noise=GREY_LEVEL;
 if (noise<0) noise=0;
 return noise;
}

/****************************************************/
/* add gaussian noise to matrix                     */
/****************************************************/
void add_gaussian_noise(float** mat,int lgth,int wdth,float var)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 if (var!=0.0) mat[i][j]=gaussian_noise(var,mat[i][j]);
}

void add(float** matr, float** mat1, float** mat2,int lgth, int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++)
   matr[i][j]=mat1[i][j]+mat2[i][j];
}

void substract(float** matr, float** mat1, float** mat2, int lgth, int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++)
   matr[i][j]=mat1[i][j]-mat2[i][j];
}

void copy(float** mat_in, float** mat_out,int lgth, int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++)
   mat_in[i][j]=mat_out[i][j];
}


void Recal_haar_step(float** image,int x0,int xf,int y0,int yf)
{
  int i,j;
  float min,max;

  min=image[x0][y0];
  max=min;
  for(i=x0;i<xf;i++)
    for(j=y0;j<yf;j++)
      {
	if(min>image[i][j]) min=image[i][j];
	if(max<image[i][j]) max=image[i][j];
      }
  for(i=x0;i<xf;i++)
    for(j=y0;j<yf;j++)
      {
	image[i][j]-=min;
	image[i][j]/=(max-min);
	image[i][j]*=255;
      }
}

void Recal_haar(float** image,int M,float** mat_tmp,int lgth,int wdth)
{
  int a,count,i,j;

  for(i=0;i<lgth;i++)
    for(j=0;j<wdth;j++)
      mat_tmp[i][j]=image[i][j];

  a=1;
  for(count=1;count<=M;count++)
    {
      Recal_haar_step(mat_tmp,lgth/(2*a),lgth/a,0,wdth/(2*a));
      Recal_haar_step(mat_tmp,lgth/(2*a),lgth/a,wdth/(2*a),wdth/a);
      Recal_haar_step(mat_tmp,0,lgth/(2*a),wdth/(2*a),wdth/a);
      a*=2;
    }
  Recal_haar_step(mat_tmp,0,lgth/(a),0,wdth/(a));
}



void haar1D(float* signal,float* work,int lgth)
{
  int i,x;

  for(i=0;i<lgth;i++)
    work[i]=signal[i];
  
  for(i=0;i<lgth/2;i++)
    {
      x=2*i;
      signal[i]=(work[x]+work[x+1])/sqrt(2);
      signal[i+lgth/2]=(work[x]-work[x+1])/sqrt(2);
    }
}

void haar2D(float** image,int lgth,int wdth)
{
  int i,j;
  float* work=fmatrix_allocate_1d(wdth);
  float* work2=fmatrix_allocate_1d(lgth);
  float* tmp2=fmatrix_allocate_1d(lgth);
  
  for(i=0;i<lgth;i++)
    haar1D(image[i],work,wdth);

  for(j=0;j<wdth;j++)
    {
      for(i=0;i<lgth;i++)
	tmp2[i]=image[i][j];
      haar1D(tmp2,work2,lgth);
      for(i=0;i<lgth;i++)
	image[i][j]=tmp2[i];
    }
  
  free_fmatrix_1d(work);
  free_fmatrix_1d(work2);
  free_fmatrix_1d(tmp2);
}

void haar2D_complete(float** image,float** haar,int M,int lgth,int wdth)
{
  int a,count;
  int i,j;
  
  for(i=0;i<lgth;i++)
    for(j=0;j<wdth;j++)
      haar[i][j]=image[i][j];
  
  a=1;
  for(count=1;count<=M;count++)
    {
      haar2D(haar,lgth/a,wdth/a);
      a*=2;
    }
}

void ihaar1D(float* signal,float* work,int lgth)
{
  int i,x;

  for(i=0;i<lgth;i++)
    work[i]=signal[i];
  
  for(i=0;i<lgth/2;i++)
    {
      x=2*i;
      signal[x]=(work[i]+work[i+lgth/2])/sqrt(2);
      signal[x+1]=(work[i]-work[i+lgth/2])/sqrt(2);
    }
}

void ihaar2D(float** image,int lgth,int wdth)
{
  int i,j;
  float* work=fmatrix_allocate_1d(wdth);
  float* work2=fmatrix_allocate_1d(lgth);
  float* tmp2=fmatrix_allocate_1d(lgth);
  
  for(i=0;i<lgth;i++)
    ihaar1D(image[i],work,wdth);

  for(j=0;j<wdth;j++)
    {
      for(i=0;i<lgth;i++)
	tmp2[i]=image[i][j];
      ihaar1D(tmp2,work2,lgth);
      for(i=0;i<lgth;i++)
	image[i][j]=tmp2[i];
    }
  
  free_fmatrix_1d(work);
  free_fmatrix_1d(work2);
  free_fmatrix_1d(tmp2);
}

void ihaar2D_complete(float** haar,float** haar_inverse,int M,
		      int lgth,int wdth)
{
  int a,count;
  int i,j;
  
  for(i=0;i<lgth;i++)
    for(j=0;j<wdth;j++)
      haar_inverse[i][j]=haar[i][j];
  
  a=1;
  for(count=1;count<=M-1;count++)
    a*=2;
  for(count=1;count<=M;count++)
    {
      ihaar2D(haar_inverse,lgth/a,wdth/a);
      a/=2;
    }
}

