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

#define TT 3000 // maximal simulation time
#define N 5000 // size of the network 

float AA[N][N],A[N][N]; // adjacency matrix
double sigma,beta,alpha,Psi[N],X; // mutation rate, infection rate, recovery rate, fitness, mean of gaussian random walk
double state[N]; // state of a individual, -1=susceptible, 1=infected, -2=recoveried
double SS,RR,I,aveI,aveR,avecount; //variables for statistics
double linum[N]; //number of infected neighbours of nodes
double count,count2,count3,judge; // number of realization, number of epidemic outbreak, times that the fitness can reach the shreshold, judgement: whether the  fitness reaches the shreshold
int tR;// timing of the start of mitigation 

// generating ER network 
void network(){
	int i,j,suml,h,hh,k;
	double tempa,tempb,st[N];;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
		    A[i][j]=0;
		    A[j][i]=0;
		}
	}
	
	
	suml=0;
 	for(i=0;i<N;i++)    
		{for(j=0;j<N;j++)
		{A[i][j]=0;}}
	while(suml<37500)
	{h=rand()%N;
 	hh=rand()%N;
 	if((h!=hh)&&(fabs(A[h][hh])<0.01))
 	{A[h][hh]=1;
 	 A[hh][h]=1;
 	 suml=suml+1;}
	}	
}

// main program
int main(){
    int i,j,jj,t,k,kk,h; // variables for loops
    double a,d,temp;// temporal variables for calculation 

    FILE* fp;
    fp=fopen("da.txt","w");

	alpha=0.1;
	X=0;
	srand(time(NULL));   
	sigma=0.03; 
	
	  	  
for(tR=0;tR<=120;tR=tR+2) 
	{

		aveI=0;
		avecount=0;
		count2=0;
		count3=0;
	for(count=0;count<100;count++) // 50 realizations for each parameters
	{
	beta=0.008;		
	judge=-1; 
	for(i=0;i<N;i++) for(j=0;j<N;j++) A[i][j]=0;
    
    network();
    for(i=1;i<N;i++)
   	  {for(j=1;j<N;j++)
  		  {AA[i][j]=A[i][j];
	    	}
		}

	
	  for(i=0;i<N;i++) // initialize
	 {       
	   linum[i]=0;	   
	 }	
	  
		 SS=0;
		 I=0;
		 RR=0;
	 	 
	 for(i=0;i<N;i++){
		  state[i]=-1;
		  Psi[i]=1;
	  }
	  for(i=0;i<N/50;i++)
	  {
	  	state[i]=1;
	  }
	  
	  for(i=0;i<N;i++)   
	  {
		if(state[i]<0)
			{
			SS++;
				}
	   	if(state[i]>0)
			{
			I++;
			}			
	  } 
	  
	 aveI=0;
	 aveR=0;
     avecount=0;
		
	 for(t=0;t<TT;t++) 
	 {	 
	 //when the mitigation starts, the infection rate is set as the half of the original
	 if(t>=tR)
	 	{
	 		beta=0.004;
		 }
		 
	 	      for(i=0;i<N;i++) // linum[i] record the number of infected neighbours of node i
	  {
	     temp=0;      
		 for (j=0;j<N;j++)  
		 {	if((AA[i][j]>0.5) && (state[j]>0))
			 {
			   temp++;			 
			 }
		 }
	     linum[i]=temp;
	  }
	  // infection process
	  for (i=0;i<N;i++)
	  {
	  	if ((state[i]<0)&&(state[i]>-1.5))
	  	{
	  		a=rand()/(1.0+RAND_MAX);
	  		if(a<1-pow((1-beta),linum[i] ))
	  		{
	  			state[i]=1;
	  			k=(int)linum[i];
	  			kk=rand() % k +1;
	  			jj=0;
	  			for(j=0;j<N;j++)
	  			{
				   if((AA[i][j]>0.5) && (state[j]>0))
	  			   {
	  			   	jj=jj+1;
	  			   	if(jj==kk)
	  			   	{
						 Psi[i]=Psi[j];
						 }
					 }
					 
	  		    }
			  }
		  }
		  else if(state[i]>=0 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<alpha/Psi[i])
		  	{
			  state[i]=-2;
			  }
			  
		  }
	  }
	 
	//statistics 
   	SS=0;
	I=0;	
	RR=0;
	temp=0;	
	  for(i=0;i<N;i++)   
	  {
		if((state[i]<0)&&(state[i]>-1.5))
			{
			SS++;
				}
		if(state[i]<-1.5)
		{
			RR++;
				}		
	   	if(state[i]>0)
			{
			I++;
			d=0;
			//random walk of fitness
			for(h=0;h<12;h++)
			{
				a=(double)((double)rand()/(double)RAND_MAX);
                d=d+a;
			}
			d=d-6;
		    d=X+d*sigma;
		    Psi[i]=Psi[i]+d;
		    if(Psi[i]<0)
		    {
		    	Psi[i]=0;
			}
			
			if(Psi[i]>20)
		    {
		    	Psi[i]=20;
			}
			
		    temp=temp+Psi[i];
			}			
	  }
	  temp=temp/I; 	

	//check if the average fitness reaches the thershold creating the second wave
	  if(temp>1.67)
	  {
	  	judge=1;
		  }	
	  
	 if(t>2500)
	 {
	 aveI=aveI+I/N;
	 aveR=aveR+RR/N;
	 avecount=avecount+1;

	 } 
	 
	 
 } 

	 // count number of outbreak
	 if(RR/N>=0.2)
	 {
	 	count2=count2+1;
	 }
	 if(judge>0)
	 {
	 	count3=count3+1;
	 }
	

}

fprintf(fp,"%d  %f  %f  %f  %f\n",tR,beta,sigma,count2/100,count3/100);	 //output

}
    return 1;
}
