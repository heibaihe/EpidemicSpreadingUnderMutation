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
double linum[N]; //probability of being infected for each node 
double linum2[N]; //number of infected neighbours of nodes
double count,count2; // number of realization, number of epidemic outbreak
double pI[N]; // probability of nodes' fitness that is transfered to the newly infected node

// generating ER network 
void network(){
	int i,j,suml,h,hh,k,kk;
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
for(sigma=0.001;sigma<=10;sigma=sigma*pow(10,0.04))	  
{	  	  
for (beta=0;beta<=0.01;beta=beta+0.0002) // loops for R0
	{

		count2=0;
	for(count=0;count<50;count++) // 50 realizations for each parameters
	{
			
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
	 	      for(i=0;i<N;i++) // linum[i] record the number of infected neighbours of node i
	  {
	     temp=0;  
		 linum[i]=1;     
		 for (j=0;j<N;j++)  
		 {	if((AA[i][j]>0.5) && (state[j]>0))
			 {
			   temp++;	
			   linum[i]=linum[i]*(1-beta*Psi[j]);		 		 
			 }
		 }
	     linum2[i]=temp;
	  }
	  // infection process
	  for (i=0;i<N;i++)
	  {
	  	if ((state[i]<0)&&(state[i]>-1.5))
	  	{
	  		a=rand()/(1.0+RAND_MAX);
	  		if(a<1-linum[i] )
	  		{
	  			//fitness of the node being infected is equal to the fitness of the node infecting it
	  			temp=0;
	  			for (j=0;j<N;j++)  
		       {	if((AA[i][j]>0.5) && (state[j]>0))
			   {
			   temp=temp+Psi[j];	
			   }
			   
		       }
		       kk=0;
		       for (j=0;j<N;j++)  
		       {	if((AA[i][j]>0.5) && (state[j]>0))
			   {
			   	if(kk==0)
			     {
				 pI[kk]=Psi[j]/temp;
				 }
				else if(kk>0)
				 {
				 	pI[kk]=Psi[j]/temp+pI[kk-1];
				 }
				 kk=kk+1; 
			   }			   
		       }
		       
		       a=rand()/(1.0+RAND_MAX);
		       if(a<pI[0])
		       {
		       	kk=1;
			   }
			   else if (a>pI[0])
			   {
			   	for(j=1;j<N;j++)
			   	{
			   		if ((a>pI[j-1])&&(a<pI[j]))
			   		{kk=j;
					   }
				   }
			   }
			   
			   
	  			state[i]=1;
	  			k=(int)linum2[i];
	  			
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
		  	if(a<alpha)
		  	{
			  state[i]=-1;
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
			
			if(Psi[i]>10)
		    {
		    	Psi[i]=10;
			}
			
		    temp=temp+Psi[i];
			}			
	  }
	  temp=temp/I; 	
	  
	 if(t>2500)
	 {
	 aveI=aveI+I/N;
	 aveR=aveR+RR/N;
	 avecount=avecount+1;

	 } 
	 
	 
 } 

	 // count number of outbreak
	 if(I/N>=0.2)
	 {
	 	count2=count2+1;
	 }
	 
	

}

fprintf(fp,"%f  %f  %f\n",beta,sigma,count2/50);	 //output

}} 
    return 1;
}
