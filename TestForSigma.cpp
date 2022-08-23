/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/
/** **************************************************************************************************** **/

#include "math.h"
//#include "iostream.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"


#define N 200 //number of strains
#define M 2000 // number of individuels
#define K 5000 // number of realizations




double F[N][2]; //fitness of a strain, F[i][0] is the intra-host fitness of strain i, F[N][1] is the inter-host fitness;
double q,qq,a,b,c,d; // temperal variables 
double sigmaPsi,sigmaPhi,sigma; // deviation of fitness
double ave1,ave2,ave3,std1,std2,aveK1,aveK2,stdK1,stdK2;// variables for statistics  
int i,j,k,ii,jj,kk,h,kkk,iii,t; // variables for loops
int strain[M]; // strian index
double p; //mutation rate
int rho; // number of iteration


int main(){

    FILE* fp;


    fp=fopen("da.txt","w");
  

    srand(time(NULL)); 

p=0.02;
rho=1000;
sigmaPhi=0.2;
//sigmaPsi=0.2;

//for(p=0;p<0.5;p=p+0.002)
//for(rho=0;rho<1000;rho=rho+5)
//for(sigmaPhi=0;sigmaPhi<0.6;sigmaPhi=sigmaPhi+0.02)
for(sigmaPsi=0;sigmaPsi<0.6;sigmaPsi=sigmaPsi+0.02)
{    
aveK1=0;
aveK2=0;
stdK1=0;
stdK2=0;

for(ii=0;ii<K;ii++)
{
	//set the fitness for strains 
	
    for(i=0;i<N;i++)
    {
    	q=0;
		qq=0; 
		for(h=0;h<12;h++)
		{
			a=(double)((double)rand()/(double)RAND_MAX);
			q=q+a;
		}
			q=q-6;
			qq=1+q*sigmaPsi;
			if(qq<0)
			{
				qq=0;
			}
			if(qq>2)
			{
				qq=2;
			 } 
			F[i][0]=qq;    	
	}
	
	for(i=0;i<N;i++)
    {
    	q=0;
		qq=0; 
		for(h=0;h<12;h++)
		{
			a=(double)((double)rand()/(double)RAND_MAX);
			q=q+a;
		}
			q=q-6;
			qq=1+q*sigmaPhi;  
			if(qq<0)
			{
				qq=0;
			}
			if(qq>2)
			{
				qq=2;
			 }
			F[i][1]=qq;
		
		
	}
	

	// initialize all individuals as a random chosen strain
	for(i=0;i<M;i++)
	{
		strain[i]=N/2;
	}

//iterations  
for(t=0;t<rho;t++)
{

	//statistics
	ave1=0;
	ave2=0;
	ave3=0;
	std1=0;
	std2=0;
	 
	for(i=0;i<M;i++)
	{
		ave1=ave1+F[strain[i]][0];
		ave2=ave2+F[strain[i]][1];
		ave3=ave3+strain[i];

	}
	ave1=ave1/M;
	ave2=ave2/M;
	ave3=ave3/M;
	
	for(i=0;i<M;i++)
	{
		std1=std1+(F[strain[i]][0]-ave1)*(F[strain[i]][0]-ave1);
		std2=std2+(F[strain[i]][1]-ave2)*(F[strain[i]][1]-ave2);

	}
	std1=sqrt(std1/M);
	std2=sqrt(std2/M);
	
	if(t==rho-1)
	{
		aveK1=aveK1+ave1/K;
		aveK2=aveK2+ave2/K;
		stdK1=stdK1+std1/K;
		stdK2=stdK2+std2/K;
		
	}
	

	
	for(i=0;i<M;i++)
	{
		a=(double)((double)rand()/(double)RAND_MAX);
		if(a<=p)
		{
			strain[i]++;
			if(strain[i]>N-1)
			{
				strain[i]=N-1;
			}
		}
		if((a>p)&&(a<=2*p))
		{
			strain[i]--;
			if(strain[i]<0)
			{
				strain[i]=0;
			}
		}	
		
	
	}
	
	//propagation of strains with higher intra-host fitness
	for(i=0;i<M;i++)
	{
		a=(double)((double)rand()/(double)RAND_MAX);
		
		if(a<0.1/F[strain[i]][0])
		{
			strain[i]=-1;
			b=-1;
			while(b<0)
			{
			j=rand() % M +1;
			if(strain[j]>=0)
			{
				strain[i]=strain[j];
				b=1;
			}
			}
		}

		
		
	}
	
	
	
}


} 

fprintf(fp,"%f	 %f	 %f	 %f	 %f  %f	 %f  %f\n",rho,p,sigmaPsi,sigmaPhi,aveK1,aveK2,stdK1,stdK2);

}
	
    return 0;
}
