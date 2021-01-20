// HMM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<cstring>
#include <sstream>
#include <algorithm>
#include <queue>
#define states 5
#define outputs 32
using namespace std;


vector<vector<long double>> A(states,vector<long double>(states,0)); //state transition matrix, n*n
vector<vector<long double>> B(states,vector<long double>(outputs)); //state observation matrix n*possible observations
vector<long double> Pi(states); //Intial state prob
vector<vector<long double>> alpha;		//T*n
vector<vector<long double>> beta;		//T*n
vector<vector<long double>> codebook;
vector<int> O;
vector<vector<int>> obs_sequences;
#define PI 3.14159265

//hamming window
vector<long double> hamming_window;

//liftering window
vector<long double> liftering_window;


long double forward(int n)
{
	long double prob_of_O_given_model=0.0;
	alpha.clear();
	//Forward_prop
	int T=O.size();
	//cout<<"\nAt alpha function T="<<T<<"\n";
	for(int i=0;i<T;i++)
	{
		alpha.push_back(vector<long double>(n,0.0));
	}
	
	for(int i=0;i<n;i++)
	{
		alpha[0][i]=(long double)(Pi[i]*B[i][O[0]]);
	}
	
	//Inductive step
	for(int t=1;t<T;t++)
	{
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				alpha[t][i]+=(long double)(alpha[t-1][j]*A[j][i]);							
			}
			alpha[t][i]=(long double)((long double)alpha[t][i]*(long double)B[i][O[t]]);
		}
	}
	
	for(int i=0;i<n;i++)
	{
		prob_of_O_given_model+=alpha[T-1][i];
	}
	return prob_of_O_given_model;
}

//Calculating Beta for backward algorithm

void Backward(int n)
{
	beta.clear();
	int T=O.size();

	for(int i=0;i<T;i++)
	{
		beta.push_back(vector<long double>(n,0.0));
	}

	for(int i=0;i<n;i++)
	{
		beta[T-1][i]=1.0;
	}

	for(int t=T-2;t>=0;t--)
	{
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				beta[t][i]+=(long double)(A[i][j]*B[j][O[t+1]]*beta[t+1][j]);
			}
		}
	}
	return;
}


//Soln to 2, Vitarbi algorithm, calculating delta

vector<vector<long double>> delta;		//T*n stores highest prob along a single path which accounts for first t observations and ending at state i given the model. delta[t+1][j]=max(delta[t][i]*Aij) * B[j][O[t+1]]
vector<int> Q; //stores state which is most optimal to occur at time t.
vector<vector<int>> sigma; // stores the a previous state which is best for current state at time t. t*n matrix
long double P_star;
void vitarbi(int n)
{
	delta.clear();
	sigma.clear();
	Q.clear();
	int T=O.size();
	//cout<<T<<"\n";
	for(int i=0;i<T;i++)
	{
		vector<long double> temp(n,0.0);
		vector<int> temp1(n,0);
		delta.push_back(temp);
		sigma.push_back(temp1);
	}
	//cout<<"Delta Matrix\n";
	for(int i=0;i<n;i++)
	{
		delta[0][i]=(long double)Pi[i]*B[i][O[0]];
		//cout<<delta[0][i]<<" ";
		sigma[0][i]=-1;
	}
	cout<<"\n";
	//recursive
	for(int t=1;t<T;t++)
	{
		for(int i=0;i<n;i++)
		{
			long double maxvalue=std::numeric_limits<long double>::min();
			int index=0;
			for(int j=0;j<n;j++)
			{
				if(delta[t-1][j]*A[j][i]>=maxvalue)
				{	//cout<<"\nEntered inside\n";
					maxvalue=(long double)(delta[t-1][j]*A[j][i]);
					index=j;
				}
			}
			delta[t][i]=(long double)(maxvalue*B[i][O[t]]);
			//cout<<delta[t][i]<<" ";
			sigma[t][i]=index;
		}
		//cout<<"\n";
	}

	long double P=std::numeric_limits<long double>::min(); //max value out of sequences ending at ith state at time T.
	int index=0;
	for(int i=0;i<n;i++)
	{
		if(delta[T-1][i]>P)
		{
			index=i;
			P=delta[T-1][i];
		}
	}
	P_star=P;
	//cout<<P_star<<"\n";
	Q.push_back(index);
	for(int t=T-2;t>=0;t--)
	{
		int index=Q[Q.size()-1];
		Q.push_back(sigma[t+1][index]);
	}
	reverse(Q.begin(),Q.end());
	return;
}

//Now we go to prob 3 which is a reestimation problem and here i need to calculate Xhi matrix and gyma matrix.
//Xhi matrix is 3D matrix T*n*n which stores prob of being in state i int time t and state j in time t+1. Xhi[t][i][j]=P[qt=Si,qt+1=Sj/O,model]
//Gyma matrix is a 2D matrix which stores prob of being in state i at time t.
//This both helps in reestimation of model.
//How?, t=1 to T sum(Xhi[t][i][j]) tells expected number of transitions from state i to j over all time period, same for gamma for t=1 to T sum(gamma[t][i]) tells expected number of time being in state i.
vector<vector<vector<long double>>> Xhi; //t*n*n matrix
vector<vector<long double>> gamma;		//t*n matrix;

void Xhi_gamma_calculations(int n)
{
	Xhi.clear();
	gamma.clear();
	int T=O.size();
	for(int i=0;i<T;i++)
	{
		Xhi.push_back(vector<vector<long double>>(n,vector<long double>(n,0.0)));
		gamma.push_back(vector<long double>(n,0.0));
	}

	vector<long double> total_t;
	for(int t=0;t<T-1;t++)
	{
		long double sum=0.0;
		for(int i=0;i<n;i++)
		{
			long double alpha_t_i=alpha[t][i];
			for(int j=0;j<n;j++)
			{
				sum+=(long double)(alpha_t_i*A[i][j]*B[j][O[t+1]]*beta[t+1][j]);
			}
		}
		total_t.push_back(sum);
	}

	//gamma calculation
	for(int t=0;t<T;t++)
	{
		long double sum1=0.0;
		for(int j=0;j<n;j++)
		{
			sum1+=(long double)(alpha[t][j]*beta[t][j]);
		}
		for(int i=0;i<n;i++)
		{
			gamma[t][i]=long double((long double)(alpha[t][i]*beta[t][i])/(long double)sum1);
		}
	}

	// Now we calculate Xhi and gamma; gamma means prob of being in state i at time t, now as we have Xhi, we can sum up Xhi[t][i][j] where j=1 to n and get gamma [t][i];
	for(int t=0;t<T-1;t++)
	{
		for(int i=0;i<n;i++)
		{
			//gamma[t][i]=0.0;
			for(int j=0;j<n;j++)
			{
				Xhi[t][i][j]=(long double)((long double)(alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j]) / (long double)total_t[t]);
				//gamma[t][i]+=Xhi[t][i][j];
			}
		}
	}
	return;
}

void reestimation(int n)
{
	//First we reestimate Pi, Expected frequency in state i at time t=1.
	int T=O.size();
	/*cout<<"Pi vector\n";
	for(int i=0;i<n;i++)
	{
		Pi[i]=gamma[0][i];
		cout<<Pi[i]<<" ";
	}*/
	cout<<"\nA Matrix\n";
	//Now reestimate A matrix;
	
	vector<long double> gamma_state;
	for(int i=0;i<n;i++)
	{
		long double sum_gamma=0.0;	//Expected no of time being in state i.
		for(int t=0;t<T;t++)
		{
			sum_gamma+=gamma[t][i];
		}
		gamma_state.push_back(sum_gamma);
	}
	
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			long double sum_xhi=0.0;		//Expected no of transitions from state i to j
			long double sum_gamma=0.0;
			for(int t=0;t<T-1;t++)
			{
				sum_xhi+=Xhi[t][i][j];
				sum_gamma+=gamma[t][i];
			}
			if(gamma_state[i]!=0)
				A[i][j]=(long double)((long double)sum_xhi/sum_gamma);
			else
				A[i][j]=0.0;
			cout<<A[i][j]<<" ";
		}
		cout<<"\n";
	}
	

	cout<<"B Matrix\n";
	//Now we reestimate B
	for(int j=0;j<n;j++)
	{
		for(int k=0;k<32;k++)
		{
			long double sum_num=0.0;
			for(int t=0;t<T;t++)
			{
				if(O[t]==k)
				{
					sum_num+=(long double)gamma[t][j];
				}
			}
			long double power=pow(10.0,-30);
			if(gamma_state[j]!=0.0)
				B[j][k]=max(power,(long double)((long double)sum_num/(long double)gamma_state[j]));
			else
				B[j][k]=power;
			cout<<B[j][k]<<" ";
		}
		cout<<"\n";
	}
	return;
}

//From frame which has max_ste we consider some frames from left and some from right of it and complete total needed range of frames.
void calculate_steady_frame(vector<long double> data,vector<long double> &steady_data, long long int middle, int range)
{
	long long int start;
	long long int end;
	if(range%2!=0)
	{
		start=middle-((range/2)+1)*80;
		end=middle+(range/2)*80;
	}
	else{
		start=middle-(range/2)*80;
		end=middle+(range/2)*80;
	}
	//start index ko check krna hain!
	if(start<0)
	{
		start=0;
		end=320+(range-1)*80;
	}
	else if(end>=data.size())
	{
		end=data.size();
		start=data.size()-(320+(range-1)*80);
	}
	if(data.size()<=7040)
	{
		start=0;
		end=data.size();
	}
	cout<<"\nStart Frame is "<<start<<" end frame is "<<end<<"\n";
	for(long long int i=start;i<end;i++)
	{
		steady_data.push_back(data[i]);
	}
	return;
}

//Vector Quantization: Here we take input as a vector of ceps coeff of a frame, find min dist of vector out of all the codebook vector and than we add obtained index in Obs sequence;

void Vector_Quantization(const vector<long double> ceps_frame)
{
	int index=0;
	long double min_error=std::numeric_limits<double>::max();
	long double tokura_weights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	for(int i=0;i<codebook.size();i++)
	{	
		long double error_sum=0;
		for(int j=0;j<12;j++)
		{
			error_sum+=(long double)tokura_weights[j]*(pow(ceps_frame[j]-codebook[i][j],2));
		}	
		if(error_sum<min_error)
		{
			min_error=error_sum;
			index=i;
		}
	}
	O.push_back(index);
}

//find Ri value for given single steady frame.
void find_Ri(vector<long double> data,long double *Ri,long long int start, long long int end,int p)
{
	//Apply hamming window on data;
	int count=0;
	for(long long int i=start;i<end;i++)
	{
		data[i]=data[i]*hamming_window[count];
		count++;
	}
	for(int i=0;i<p+1;i++)
	{
		for(long long int j=start;j<end-i;j++)
		{
			Ri[i]+=data[j]*data[j+i];
		}
		//cout<<"Frame from "<<start<<"  "<<end<<"  "<<Ri[i]<<"\n";
	}
	//cout<<"\n";
	return;
}


//Does LPC analysis on given data and finds ais and cis subsequently for each frame.
void LPC_analysis(vector<long double> data,int p)
{
	long long int window_size=320;
	long long int end=window_size;
	long long int sliding_window=80;
	long long int start=0;
	int frame_count=0;
	while(end<data.size())
	{
		long double *R= new long double[p+1];
		for(int i=0;i<=p;i++)
			R[i]=0;
		find_Ri(data,R,start,end,p);
		long double *E=new long double [p+1];
		long double *ai = new long double[p+1];
		long double *ci = new long double[p+1];
		E[0]=R[0];
		long double *K= new long double[p+1];
		long double** alpha = new long double*[p+1];
		for(int i = 0; i < p+1; ++i)
			alpha[i] = new long double[p+1];
		for(int i=1;i<=p;i++)
		{
			long double sum_k=0;
			for(int j=1;j<=i-1;j++)
			{
				sum_k+=alpha[i-1][j]*R[i-j];
			}
			K[i]=(R[i]-sum_k)/E[i-1];
			alpha[i][i]=K[i];
			for(int j=1;j<=i-1;j++)
			{
				alpha[i][j]=alpha[i-1][j] - K[i]*alpha[i-1][i-j];
			}
			E[i]=(1-K[i]*K[i])*E[i-1];
		}
		for(int j=1;j<=p;j++)
		{
			ai[j]=alpha[p][j];
			//cout<<ai[j]<<" ";
		}
		//cout<<"\n";
		
		//counting ci values;
		ci[0]=log(R[0]*R[0]);
		for(int i=1;i<=p;i++)
		{
			long double temp_sum=0;
			for(int j=1;j<i;j++)
			{
				long double j_i = (long double)j/(long double)i;
				temp_sum += (ci[j]*ai[i-j])*j_i;
			}
			temp_sum+=ai[i];
			ci[i]=temp_sum;
			//cout<<ci[i]<<" ";
		}
		vector<long double> cepstral_coeffecient_for_frame;
		for(int j=1;j<=p;j++)
		{
			ci[j]=ci[j]*liftering_window[j];
			//cout<<ci[j]<<"  ";
			cepstral_coeffecient_for_frame.push_back(ci[j]);
		}
		//cout<<"\n";
		Vector_Quantization(cepstral_coeffecient_for_frame);
		frame_count++;
		start+=sliding_window;
		end+=sliding_window;
		delete[] E;
		delete[] R;
		delete[] ci;
		delete[] ai;
		delete[] K;
		for(int i = 0; i < p+1; ++i)
			delete[] alpha[i];
		delete[] alpha;

	}
	cout<<"Size of observation sequence: "<<frame_count<<"\n";
}

void intialize_model(int n,int m,int i)
{
	long double Prob=(i==0)?0.8:0.9;
	//first we intialize A matrix with bakers model 0.8,0.2 
	for(int i=0;i<n-1;i++)
	{
		A[i][i]=Prob;
		A[i][i+1]=(long double)(1.0-Prob);
	}
	A[n-1][n-1]=1;
	//Now we intialize B matrix where each event has equal probability to occur.
	long double prob=(long double)(1.0/m);
	cout<<"Intialization B matrix\n";
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			B[i][j]=prob;
			cout<<B[i][j]<<" ";
		}
		cout<<"\n";
	}

	//Now we intialize Pi Matrix where starting has 1 prob and rest 0.
	Pi[0]=1.0;
	for(int i=1;i<n;i++)
	{
		Pi[i]=0.0;
	}
	return;
}

//Find best models from all the models and avg them out. Best as per max P* values!

void retrive_best_models_for_a_digit(int best)
{
	//Intialize A,B matrix to 0.
	for(int i=0;i<states;i++)
	{
		for(int j=0;j<states;j++)
		{
			A[i][j]=0.0;
		}
	}

	for(int i=0;i<states;i++)
	{
		for(int j=0;j<outputs;j++)
		{
			B[i][j]=0.0;
		}
	}

	//I need to collect data from P_star file and store them in a array and then pick best n values from it and corressponding models from models file.
	ifstream model_dump;
	ifstream Pstar_dump;
	Pstar_dump.open("pstar_dump_file.txt");

	//store values in vector. we can make it very effecient by using min heap of size K.
	vector<long double> p_star;
	priority_queue<long double> pq;
	long double temp;
	while(Pstar_dump>>temp)
	{
		pq.push(temp);
	}

	cout<<"\nPstar of best "<<best<<" models: \n";
	for(int i=0;i<best;i++)
	{
		cout<<pq.top()<<"\t";
		p_star.push_back(pq.top());
		pq.pop();
	}
	cout<<"\n";
	model_dump.open("model_dump_file.txt");
	int flag=0;
	int counter=0;
	while(model_dump>>temp)
	{
		for(int k=0;k<p_star.size();k++)
		{
			//cout<<temp<<"\t";
			if(temp==p_star[k])
			{
				flag=1;
				counter++;
				cout<<"Inner side"<<temp<<"\t";
				for(int i=0;i<states;i++)
				{
					for(int j=0;j<states;j++)
					{
						model_dump>>temp;
						A[i][j]+=temp;
					}
				}

				for(int i=0;i<states;i++)
				{
					for(int j=0;j<outputs;j++)
					{
						model_dump>>temp;
						B[i][j]+=temp;
					}				
				}
				break;
			}
		}
		if(counter==best)
			break;
		if(!flag)
		{
				for(int i=0;i<states;i++)
				{
					for(int j=0;j<states;j++)
					{
						model_dump>>temp;
					}
				}

				for(int i=0;i<states;i++)
				{
					for(int j=0;j<outputs;j++)
					{
						model_dump>>temp;
					}				
				}
		}
		flag=0;
	}
	cout<<"\n";
	ofstream final_model_file;
	final_model_file.open("final_model_file.txt",std::ios_base::app);
	long double mult=(long double)((long double)1.0/(long double)best);
	//Printing A matrix after averaging out.
	for(int i=0;i<states;i++)
	{
		for(int j=0;j<states;j++)
		{
			A[i][j]*=mult;
			cout<<A[i][j]<<" ";
			if(best==19)
				final_model_file<<A[i][j]<<"\t";
		}
		if(best==19)
			final_model_file<<"\n";
		cout<<"\n";
	}
	//long double power=pow(10.0,-30);
	//Printing B matrix
	for(int i=0;i<states;i++)
	{
		for(int j=0;j<outputs;j++)
		{
			B[i][j]*=mult;
			cout<<B[i][j]<<" ";
			if(best==19)
				final_model_file<<B[i][j]<<"\t";
		}
		if(best==19)
			final_model_file<<"\n";
		cout<<"\n";
	}
	Pstar_dump.close();
	model_dump.close();
	final_model_file.close();


	if( remove( "model_dump_file.txt" ) != 0 )
		perror( "Error deleting file" );
	if( remove("pstar_dump_file.txt") != 0)
		perror( "Error deleting file" );
	return;
}

//Stores model and Pstar for all files to get best value from all files in digit.
void store_pre_models(int n)
{
	ofstream model_dump;
	model_dump.open("model_dump_file.txt",std::ios_base::app);
	ofstream Pstar_dump;
	Pstar_dump.open("pstar_dump_file.txt",std::ios_base::app);
	if(!Pstar_dump)
	{ 
       cout<<"Error in opening Pstar store file!!!"; 
       return; 
	}
	if(!model_dump)
	{ 
       cout<<"Error in opening Model store file!!!"; 
       return; 
	}
	Pstar_dump<<P_star<<"\n";
	model_dump<<P_star<<"\n";
	//Store model in model file.
	//Store A matrix.
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			model_dump<<A[i][j]<<"\t";
		}
		model_dump<<"\n";
	}
	model_dump<<"\n";
	//Store B matrix
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<outputs;j++)
		{
			model_dump<<B[i][j]<<"\t";
		}
		model_dump<<"\n";
	}
	model_dump<<"\n";
	//We don't need to store Pi Matrix as it never changes.
	model_dump.close();
	Pstar_dump.close();
}

//Stores data for all training files in intial training to look out how things work.
void store_data(string s,int n)
{
	ofstream model_dump;
	model_dump.open("data_dump_file_0.txt",std::ios_base::app);
	if(!model_dump)
	{ 
       cout<<"Error in opening Model store file!!!"; 
       return; 
	}
	//Store P* value
	model_dump<<s;
	model_dump<<"\n";
	for(int i=0;i<O.size();i++)
	{
		model_dump<<O[i]<<" ";
	}
	model_dump<<"\n";
	model_dump<<P_star<<"\n";
	//Store State sequence
	for(int i=0;i<Q.size();i++)
	{
		model_dump<<Q[i]<<" ";
	}
	model_dump<<"\n";
	//Store Pi matrix,
	for(int i=0;i<n;i++)
	{
		model_dump<<Pi[i]<<" ";
	}
	model_dump<<"\n";
	//Store A matrix,
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			model_dump<<A[i][j]<<" ";
		}
		model_dump<<"\n";
	}
	//Store B matrix.
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<outputs;j++)
		{
			model_dump<<B[i][j]<<" ";
		}
		model_dump<<"\n";
	}
	model_dump.close();
	return;
}

//funtions counts STE and ZCR per each frame and stores in corressponding vectors, return no of frames, it takes O(n) time complexity. Uses sliding window.

long long int Count_STE_ZCR(const vector<long double> data,vector<long double> &ste_per_frame)
{
	long double sum=0,zcr=0;											
	vector<long double> prefix_sum_energy;							//vector to store prefix sum of "square of amplitude" at each point in complete wave.
	long long int max_ste_index=0;
	  
	for(unsigned long long int i=0;i<data.size();i++)
    {
																//itereate through every data point and counts zcr and energy and sums with prefix of it
																// and stores in corresponding vectors.
        sum+=data[i]*data[i];
		prefix_sum_energy.push_back(sum);
    }
	
	//Now, counts STE and ZCR for each frame.
	int window_size=320,slide_factor=80;										//window of each frame and how much it should slide for next frame.
	long long int no_frames = (data.size()-window_size)/slide_factor + 1;		// counting total number of frames present in data.
	long long int start=0;
	long double max_ste_till_now=(long double)(prefix_sum_energy[start+window_size-1]/window_size);
	ste_per_frame.push_back((long double)(prefix_sum_energy[start+window_size-1]/window_size));
	start+=slide_factor;
	for(unsigned long long int i=1;i<no_frames;i++)
	{
		long double STE_frame=(long double)(prefix_sum_energy[start+window_size-1]-prefix_sum_energy[start-1])/window_size;
		ste_per_frame.push_back(STE_frame);		//storing STE and ZCR values for each frame
		if(max_ste_till_now<STE_frame)			//by gradually sliding by sliding factor.
		{
			max_ste_till_now=STE_frame;
			max_ste_index=i;
		}
		start+=slide_factor;
	}
	return max_ste_index;
}

void train_model(int train_time)
{
			for(int epoch=0;epoch<train_time;epoch++)
			{
				cout<<"Epoch: "<<epoch<<"\n";
				vitarbi(states);
				cout<<"P* value: "<<P_star<<"\n";
				cout<<"State Sequence: \n";
				for(int i=0;i<Q.size();i++)
				{
					cout<<Q[i]<<" ";
				}
				cout<<"\n";
				//now i need to calcuate alpha and beta.
				forward(states);
				Backward(states);
				Xhi_gamma_calculations(states);
				reestimation(states);
			}
}

int _tmain(int argc, _TCHAR* argv[])
{
	//I need to read codebook and store it in 2 dimensional vector. Codebook size is 32*12
	
	//Intializing Codebook;
	fstream codebook_file;
	codebook_file.open("codebook.txt");
	int counterr=0;
	vector<long double> codebook_single;
	//cout<<"Cepstral coeffecients from refrence files are: \n\n";
	long double tempp;
	while(codebook_file>>tempp)
	{
		codebook_single.push_back(tempp);
		//cout<<tempp<<"  ";
		counterr++;
		if(counterr%12==0)
		{
			//cout<<"\n";
			codebook.push_back(codebook_single);
			codebook_single.clear();
			counterr=0;
		}
	}
	codebook_file.close();

	int no_cepstral_coeffecients=12;

	//intializing hamming window
	for(int i=0;i<320;i++)
	{
		long double temp =(long double)( 0.54-(0.46*cos((2*PI*i)/319)));
		hamming_window.push_back(temp);
	}

	//initializing liftering window
	liftering_window.push_back(-1);
	for(int i=1;i<=no_cepstral_coeffecients;i++)
	{
		long double temp=(long double)(1+(no_cepstral_coeffecients/2)*sin((PI*i)/no_cepstral_coeffecients));
		liftering_window.push_back(temp);
	}
	/*for(int digit=0;digit<=9;digit++)
	{
		string s1="194101013_";
		stringstream ss;
		ss << digit;
		string str = ss.str();
		s1+=str;
		s1+="_";*/
	//we create model for each recording file.
	for(int f=1;f<=20;f++)
	{
		string s1="204101024_0_";
		//string s2=s1;
		stringstream ss;
		ss << f;
		string str = ss.str();
		s1+=str;
		string s3=".txt";
		s1+=s3;
		cout<<s1<<"\n";
		//Read the input file
		fstream data_file;
		data_file.open(s1);
		//I need to apply DC_Shift, Normalization and then start and end marker, after it i need to apply LPC and find Cepstrals for each frame than i need to do VQ and get Obs sequence.

		//DC_shift calculation considering initial of silence in file.
		long double temp;
		long long int counter_data = 12000;		//total sample size for 1 sec of silence in consideration.
		long long int temp_counter=counter_data;
		long double dc_shift=0.0;
		//long double silence_STE=0;
		//long double max_silence_STE=INT_MIN;
		int frame_size_counter=320;
		while(temp_counter--)
		{
			data_file>>temp;
			dc_shift+=temp;
		}

		dc_shift=dc_shift/counter_data;

		//applying DC_shift and Normalization to whole data!
		vector<long double> data;
		long double amp_req=10000;
		long double max_amp=-1*std::numeric_limits<long double>::max();
		while(data_file>>temp)
		{

			if(abs(temp-dc_shift)>max_amp)
			{
				max_amp=abs(temp-dc_shift);
			}
			data.push_back(temp-dc_shift);
		}
		data_file.close();
		cout<<"Data size: "<<data.size()<<"\n";
		
		//Normalization
		long double max_amplitude_needed=10000;
		long double ratio_amplitude=max_amplitude_needed/max_amp;
		for(unsigned long long int i=0;i<data.size();i++)
		{
			data[i]=data[i]*ratio_amplitude;
		}

		vector<long double> ste_per_frame;
		
		//Finding STE for sliding frames.
		long long int max_STE_index=Count_STE_ZCR(data,ste_per_frame);
		max_STE_index*=80;
		cout<<"Max_STE_Index"<<max_STE_index<<"\n";
		int range_data=90;
		vector<long double> steady_data;
		calculate_steady_frame(data,steady_data,max_STE_index,range_data);
		cout<<"Steady data size: "<<steady_data.size()<<"\n";
		O.clear();
		LPC_analysis(steady_data,no_cepstral_coeffecients);
		
		cout<<"Observation Sequence: \n";
		for(int i=0;i<O.size();i++)
		{
			cout<<O[i]<<" ";
		}
		cout<<"\nobservation sequence size: "<<O.size()<<"\n";
		
		obs_sequences.push_back(O);
		
		char v;
		//cin>>v;
		
		for(int i=0;i<2;i++)
		{
			intialize_model(states,outputs,i);
			train_model(40);
			//store in file model and P values.
			store_pre_models(states);
			//store_data(s1,states);
		}
	}
	//After storing all P* values i need to find best 4 models out of them and avg them out.
	retrive_best_models_for_a_digit(20);
	///*ofstream pre_model;
	//ofstream pre_pstar;
	//pre_model.open("model_dump_file.txt");
	//pre_pstar.open("pstar_dump_file.txt");
	//pre_model.close();
	//pre_pstar.close();*/
	vector<vector<long double>> pre_A;
	vector<vector<long double>> pre_B;
	for(int i=0;i<states;i++)
	{
		vector<long double> temp;
		for(int j=0;j<states;j++)
		{
			temp.push_back(A[i][j]);
			//cout<<A[i][j]<<" ";
		}
		pre_A.push_back(temp);
		temp.clear();
		//cout<<"\n";
	}
	
	for(int i=0;i<states;i++)
	{
		vector<long double> temp;
		for(int j=0;j<outputs;j++)
		{
			temp.push_back(B[i][j]);
			//cout<<B[i][j]<<" ";
		}
		pre_B.push_back(temp);
		temp.clear();
		//cout<<"\n";
	}
	//cout<<"\n";
	//Now after i find best 4 and take that as inital model and again do training for all the obs sequences.
	for(int k=0;k<obs_sequences.size();k++)
	{
		for(int i=0;i<states;i++)
		{
			for(int j=0;j<states;j++)
			{
				A[i][j]=pre_A[i][j];
				//cout<<A[i][j]<<" ";
			}
			//cout<<"\n";
		}

		for(int i=0;i<states;i++)
		{
			for(int j=0;j<outputs;j++)
			{
				B[i][j]=pre_B[i][j];
				//cout<<B[i][j]<<" ";
			}
			//cout<<"\n";
		}
		O.clear();
		O=obs_sequences[k];
		cout<<"Obs Sequence: "<<k<<" "<<O.size()<<"\n";
		for(int m=0;m<O.size();m++)
		{
			cout<<O[m]<<" ";
		}
		cout<<"\n";
		char a;
		cin>>a;
		train_model(30);
		if(P_star!=2.22507e-308 || P_star!=0)
			store_pre_models(states);
	}
	retrive_best_models_for_a_digit(19);
	vitarbi(states);
	cout<<"P* value: "<<P_star<<"\n";
	cout<<"State Sequence: \n";
	for(int i=0;i<Q.size();i++)
	{
		cout<<Q[i]<<" ";
	}
	return 0;
}

