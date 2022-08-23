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


#define TT 40000 // maximal simulation time
#define N 5000 // size of the network 



int state[N]; // state of a individual, 0=susceptible, 1=exposed, 2=miled symptoms, 3=severe symptoms, 4=critical symptoms, 5=hospitalized, 6=ventilated, 7=recovered, 8=deceased, 9=asymptomatic, 10=pre-symptomatic, 11=vaccinated
int state2[N]; // =1 if the exposed individual will become asymptomatic, =2 if the individual will become symptomatic
int state3[N]; // =-1 if the individual has been infected less than 2 times, =1 if the individual has been infected twice
double linum[N],linum2[N]; // number of neighbors those can infect a node, and the probability of a node that not being infected by its neighbors 
double linum3[N]; // the probability of a node that not being reinfected by its neighbors 
int linum4[N][50],infe[N]; // index of neighbors those can reinfect a node, number of neighbors those can reinfect a node
float A[N][N]; // adjacency matrix
double beta,sigma; // infection rate, mutation rate
double pm; // parameter controlling the reinfection rate according to the fitness
double Psi[N],Psi2[N]; // fitness of a infected node, ability of reinfection of a variant according to fitness
double count,count2; //number realizations, counts of the reemergence
double SS,EE,IM,IS,IC,H,V,R,D,AS,PS,I,Re; //variables for statistics

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
    int i,j,jj,k,t,tt,kk,kkk,h; // variables for loops
    double a,b,c,d,temp,temp2; //temporal variables for calculation 

    FILE* fp;

    fp=fopen("da.txt","w");

	
	pm=0.3;

	srand(time(NULL));
	for(i=0;i<N;i++) for(j=0;j<N;j++) A[i][j]=0;
    network();

for(sigma=0;sigma<=0.25;sigma=sigma+0.1) // loops for sigma
{

for (beta=0.0;beta<0.106;beta=beta+0.005) // loops for R0
	{


	count2=0;
	for(count=0;count<50;count++) // 50 realizations for each parameters
	{
			
	for(i=0;i<N;i++) for(j=0;j<N;j++) A[i][j]=0;
    
    network();      
	
	  for(i=0;i<N;i++) // initialize
	 {       
	   linum[i]=0;	   
	 }	
	 
	 
		 SS=0;
		 EE=0;
		 I=0;
		 R=0;
		 Re=0;
	 
	 
	 for(i=0;i<N;i++){
		  state[i]=0;
		  state2[i]=0;
		  state3[i]=-1;
		  Psi[i]=0;
		  Psi2[i]=pow(Psi[i],10)/(pow(pm,10)+pow(Psi[i],10));
		  infe[i]=1;
	  }
	  for(i=0;i<N/500;i++)
	  {
	  	state[i]=1;
	  	Psi[i]=0;
	  	Psi2[i]=pow(Psi[i],10)/(pow(pm,10)+pow(Psi[i],10));
	  	
	  	a=rand()/(1.0+RAND_MAX);
	  		if(a<0.3)
			{
 				state2[i]=1;
			}
			  else 
			  {
			  	state2[i]=2;
			  }
	  }
	  
	  for(i=0;i<N;i++)   
	  {
		if(state[i]==0)
			{
			SS++;
				}
	   	if(state[i]==1)
			{
			EE++;
			}			
	  } 
	  


		
	 for(t=0;t<TT;t++) 
	 {
	 
	 	      for(i=0;i<N;i++) // linum[i] record the number of infected neighbours of node i
	  {
	     temp=0;  
		 temp2=0;  
		 linum2[i]=1;
		 linum3[i]=1;  
		 for(j=0;j<50;j++)
		 {
		 	linum4[i][j]=-1;
		 }
		 jj=0;
		 for (j=0;j<N;j++)  
		 {	if((A[i][j]>0.5) && ((state[j]==9)||(state[j]==10)||(state[j]==12)))
			 {
			  linum2[i]=linum2[i]*(1-beta/96);
			   temp++;	
			   a=rand()/(1.0+RAND_MAX);
			   if(a<Psi2[j])
			   {
			   	temp2++;
				//linum3[i]=linum3[i]*(1-beta*Reco[j]/96);	
				linum3[i]=linum3[i]*(1-beta/96);
				linum4[i][jj]=j;	
				jj++;
					}		 
			 }
		 }
	     linum[i]=temp;
	     infe[i]=temp2;
	  }
	  
	  // dynamics of the infection cycle
	  for (i=0;i<N;i++)
	  {
	  	if (state[i]==0)
	  	{
	  		a=rand()/(1.0+RAND_MAX);
	  	if(a<1-linum2[i])
	  		{
	  			state[i]=1;
	  			
	  			a=rand()/(1.0+RAND_MAX);
	  			if(a<0.3)
	  			{
	  				state2[i]=1;
				}
				  else 
				  {
				  	state2[i]=2;
				  }
				  					 	 			   
	  			kk=(int)linum[i];
	  			kkk=rand() % kk +1;
	  			jj=0;
	  			for(j=0;j<N;j++)
	  			{
				   if((A[i][j]>0.5) && ((state[j]==9)||(state[j]==10)||(state[j]==12)))
	  			   {
	  			   	jj=jj+1;
	  			   	if(jj==kkk)
	  			   	{
						 Psi[i]=Psi[j];
						 Psi2[i]=pow(Psi[i],10)/(pow(pm,10)+pow(Psi[i],10));
						 }
					 }
					 
	  		    }	  		    
			  }
		  }
		  else if(state[i]==1 )
		  {
		  	if(state2[i]==1)
		  	{
			  a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(4.0*96))
		  	{
			  state[i]=9;
			  }
		    }
		    
		    if(state2[i]==2)
		  	{
			  a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(3.0*96))
		  	{
			  state[i]=10;
			  }
		    }
		  }
		  else if(state[i]==9 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(6.0*96))
		  	{
			  state[i]=7;
			  }
		  }
		  else if(state[i]==10 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(2.0*96))
		  	{
		  		b=rand()/(1.0+RAND_MAX);
		  		b=b;
		  		if(b<55/(70.0))
		  		{
				  state[i]=2;
				  c=rand()/(1.0+RAND_MAX);
				  if(c<-1)
				  {
				  	state[i]=12;
				  }
				  }
				if((b<65/(70.0))&&(b>55/(70.0)))
		  		{
				  state[i]=3;
				  }
				if((b>65/(70.0)))
		  		{
				  state[i]=4;
				  }
			  
			  }
		  }
		  else if(state[i]==2 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(5.0*96))
		  	{
			  state[i]=7;
			  }
		  }
		  else if(state[i]==12 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(5.0*96))
		  	{
			  state[i]=7;
			  }
		  }
		  else if(state[i]==3 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(4.0*96))
		  	{
			  state[i]=5;
			  }
		  }
		  else if(state[i]==4 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(3.0*96))
		  	{
			  state[i]=6;
			  }
		  }
		  else if(state[i]==5 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(11.0*96))
		  	{
		  		b=rand()/(1.0+RAND_MAX);
		  		if(b<0.85)
		  		{
				  state[i]=7;
				  }
			    if(b>0.85)
		  		{
				  state[i]=8;
				  }
			  }
		  }
		  else if(state[i]==6 )
		  {
		  	a=rand()/(1.0+RAND_MAX);
		  	if(a<1/(13.0*96))
		  	{
		  		b=rand()/(1.0+RAND_MAX);
		  		if(b<0.5)
		  		{
				  state[i]=7;
				  }
			    if(b>0.5)
		  		{
				  state[i]=8;
				  }
			  }
		  }	
		  else if((state[i]==7)||(state[i]==11))
		  {
		  	a=rand()/(1.0+RAND_MAX);
	  	if((a<1-linum3[i])&&(state3[i]<0))
	  		{
	  			state[i]=1;
	  			state3[i]=1;
	  			a=rand()/(1.0+RAND_MAX);
	  			if(a<0.3)
	  			{
	  				state2[i]=1;
				}
				  else 
				  {
				  	state2[i]=2;
				  }
		   
	  			kk=(int)infe[i];
	  			kkk=rand() % kk ;
	  			jj=0;
				Psi[i]=Psi[linum4[i][kkk]];
				Psi2[i]=pow(Psi[i],10)/(pow(pm,10)+pow(Psi[i],10));	  		    
			  }			  
			}	  
	  }
	  
	//statistics 	 
     SS=0;
	 I=0;	
	 EE=0;
	 IM=0;
	 IS=0;
	 IC=0;
	 H=0;
	 V=0;
	 R=0;
	 D=0;

	 temp=0;
	  for(i=0;i<N;i++)   
	  {
		if(state[i]==0)
			{
			SS++;
				}
		if((state[i]==2)||(state[i]==3)||(state[i]==4)||(state[i]==5)||(state[i]==6)||(state[i]==12))
		{
			I++;
		}
			if(state[i]==1)
			{
			EE++;
				}
			if(state[i]==7)
			{
			R++;
				}
			if(state[i]==8)
			{
			D++;
				}
			if(state[i]==5)
			{
			H++;
				}
			if(state[i]==6)
			{
			V++;
				}
							
	   	if((state[i]==1)||(state[i]==2)||(state[i]==3)||(state[i]==4)||(state[i]==5)||(state[i]==6)||(state[i]==9)||(state[i]==10)||(state[i]==12))
			{
			d=0;
			//random walk of fitness
			for(h=0;h<12;h++)
			{
				a=(double)((double)rand()/(double)RAND_MAX);
                d=d+a;
			}
			d=d-6;
		    d=0+d*sigma;
		    Psi[i]=Psi[i]+d/96;
		    
		    if(Psi[i]<0)
		    {
		    	Psi[i]=0;
			}
			
			if(Psi[i]>1)
		    {
		    	Psi[i]=1;
			}
			Psi2[i]=pow(Psi[i],10)/(pow(pm,10)+pow(Psi[i],10));
			
		    temp=temp+Psi2[i];

	    
			}			
	  }
	

	 if((I<0.5)&&(EE<0.5))
	 {t=t+40000;
	 }
	 
	 } 
	 
	 for(i=0;i<N;i++)   
	 {
	 	if(state3[i]>0)
		 {
		 	Re=Re+1;
		 }
	 }
	 
	 // count number of reemergence
	 if(Re>500)
	 {
	 	count2=count2+1;
	 }
	 


}

fprintf(fp,"%f  %f  %f\n",beta,sigma,count2/50);	//output


	  
}}
    return 1;
}
