// HMM_testing.cpp : Defines the entry point for the console application.
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

vector<vector<long double>> A(states,vector<long double>(states,0.0)); //state transition matrix, n*n
vector<vector<long double>> B(states,vector<long double>(outputs,0.0)); //state observation matrix n*possible observations
vector<long double> Pi(states); //Intial state prob
vector<vector<long double>> alpha;		//T*n
vector<vector<long double>> codebook;
vector<int> O;
#define PI 3.14159265

//hamming window
vector<long double> hamming_window;

//liftering window
vector<long double> liftering_window;

//Calculating alpha values.
long double forward(int n)
{
	long double prob_of_O_given_model=0.0;
	alpha.clear();
	//Forward_prop
	int T=O.size();
	//cout<<"\n"<<T<<"\n";
	for(int i=0;i<T;i++)
	{
		alpha.push_back(vector<long double>(n,0.0));
	}
	
	for(int i=0;i<n;i++)
	{
		alpha[0][i]=(long double)(Pi[i]*B[i][O[0]]);
	}
	
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

vector<vector<long double>> all_A;
vector<vector<long double>> all_B;

void store_all_models()
{
	ifstream data_model;
	data_model.open("final_model_m_data_40_19avg.txt");
	long double temp;
	for(int m=0;m<10;m++)
	{
		for(int i=0;i<states;i++)
		{
			vector<long double> temp_vect;
			for(int j=0;j<states;j++)
			{
				data_model>>temp;
				temp_vect.push_back(temp);
			}
			all_A.push_back(temp_vect);
			temp_vect.clear();
		}

		for(int i=0;i<states;i++)
		{
			vector<long double> temp_vect;
			for(int j=0;j<outputs;j++)
			{
				data_model>>temp;
				temp_vect.push_back(temp);
			}
			all_B.push_back(temp_vect);
			temp_vect.clear();
		}
	}
	data_model.close();
	/*for(int m=0;m<10;m++)
	{
		int start=states*m;
		for(int i=start;i<start+states;i++)
		{
			for(int j=0;j<states;j++)
			{
				cout<<all_A[i][j]<<"\t";
			}
			cout<<"\n";
		}

		for(int i=start;i<start+states;i++)
		{
			for(int j=0;j<outputs;j++)
			{
				cout<<all_B[i][j]<<"\t";
			}
			cout<<"\n";
		}
	}*/
	return;
}

//Intialize model for each digit
void initialize_model(int digit)
{
	//cout<<"Intializing model for digit "<<digit<<"\n";
	int start=(int)(digit*states);
	int counter=0;
	for(int i=start;i<start+states;i++)
	{
		for(int j=0;j<states;j++)
		{
			A[counter][j]=all_A[i][j];
			//cout<<A[counter][j]<<"\t";
		}
		//cout<<"\n";
		counter++;
	}
	counter=0;
	for(int i=start;i<start+states;i++)
	{
		for(int j=0;j<outputs;j++)
		{
			B[counter][j]=all_B[i][j];
			//cout<<B[counter][j]<<" ";
		}
		//cout<<"\n";
		counter++;
	}
	Pi[0]=1.0;
	for(int i=1;i<states;i++)
	{
		Pi[i]=0.0;
	}
	return;
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

	store_all_models();

	//Inital prerecorded file testing.
	for(int f=21;f<=30;f++)
	{
		//string s1="204101024_m_11.txt";
		string s1="204101024_9_";	
		stringstream ss;
		ss << f;
		string str = ss.str();
		s1+=str;
		string s3=".txt";
		s1+=s3;
		cout<<s1<<"\n";
		//Read the input file
		fstream data_file;
		//system("Recording_Module.exe 3 204101024_m_11.wav 204101024_m_11.txt");
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
		long double max_prob=std::numeric_limits<long double>::min();
		int ans=0;
		for(int m=0;m<10;m++)
		{
			initialize_model(m);
			long double temp_prob=forward(states);
			cout<<"Probability for digit: "<<m<<" is: "<<temp_prob<<"\n";
			if(temp_prob>max_prob)
			{
				max_prob=temp_prob;
				ans=m;
			}
		}
		cout<<" Digit Predicted is: "<<ans<<"\n";
	}

	return 0;
}

