/*
  This is the program for k-means algorithm.
  Set 3 parameters: n, d, k.
  Reona SEKINO
*/

#include<iostream>
#include<random>
#include<math.h>
#include<fstream>

using namespace std;

long double centroid(int degree, long double *trgt, long double *crtn);//calc centroid

int main(int argc, char **argv){

  if(argc != 4){//parameter error
    cout << "Assign parameter correctly" << endl;
  }
  
  int N = atoi(argv[1]);//data number
  int n = 0;//data count
  int D = atoi(argv[2]);//degree number
  int d = 0;//degree count
  int K = atoi(argv[3]);//classter number
  int k = 0;//classter count

  long double x[N][D];//data set

  //generate random data set
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> randf(0, 1000000);
  n = 0;
  while(n < N){
    d = 0;
    while(d < D){
      long double tmp = 0.0;
      tmp = randf(mt);
      x[n][d] = tmp / 1000000;
      //cout << x[n][d] << endl;//debug
      d++;
    }
    n++;
  }
  //end
  
  long double T[K][D];//represent nodes
  long double T_tmp[K][D];//tmp
  int S[N]; //class label

  clock_t start = clock();//start time
  
  //select k nodes
  int tpre = 0;
  tpre = N / K;
  k = 0;
  while(k < K){
    d = 0;
    while(d < D){
      T[k][d] = x[tpre * k][d];
      //cout << k << " " <<  d  << " " << T[k][d] << endl;//debug
      d++;
    }
    k++;
  }
  //end

  int count = 0;//repeat count
  long double mse = 0;//mean squared error
  int flag = 1;//termination flag
  while(flag){
    int S_cnt[K];//factor number of each class count
    k = 0;
    while(k < K){
      S_cnt[k] = 0;//initialize
      k++;
    }
  
    //reclasster all nodes
    n = 0;
    while(n < N){
      S[n] = 0;
      long double min = centroid(D, x[n], T[0]);
      k = 1;
      while(k < K){
      	long double tmp = centroid(D, x[n], T[k]);
	if(tmp < min){
	  min = tmp;
	  S[n] = k;
	  //cout << n << " " << S[n] << endl;//debug
	}
	k++;
      }
      S_cnt[S[n]] = S_cnt[S[n]] + 1;
      n++;
    }
    //end

    //copy T to T_tmp and initialize T
    k = 0;
    while(k < K){
      d = 0;
      while(d < D){
	T_tmp[k][d] = T[k][d];
	T[k][d] = 0;
	//cout << T_tmp[k][d] << " " << T[k][d] << endl;//debug
	d++;
      }
      k++;
    }
    //end
  
    //reset represent nodes
    n = 0;
    while(n < N){
      d = 0;
      while(d < D){
	T[S[n]][d] += x[n][d];
	d++;
      }
      n++;
    }
    k = 0;
    while(k < K){
      d = 0;
      while(d < D){
	T[k][d] = T[k][d] / S_cnt[k];
	//cout << T[k][d] << endl;//debug
	d++;
      }
      k++;
    }
    //end

    //calc mse
    mse = 0;
    n = 0;
    while(n < N){
      mse += pow(centroid(D, x[n], T[S[n]]),2);
      n++;
    }
    mse = mse / N;
    //end

    //termination judgement
    flag = 0;
    k = 0;
    while(k < K){
      d = 0;
      while(d < D){
	if(T_tmp[k][d] != T[k][d]){
	  flag = 1;
	}
	d++;
      }
      k++;
    }
    //end

    count++;//count increase
  }

  clock_t end = clock();//end time

  //report
  cout << "--- Result --- " << endl << endl;
  cout << "Data Ammount     : " << N << endl;
  cout << "Degree Ammount   : " << D << endl;
  cout << "Classter Ammount : " << K << endl;
  cout << "Time(microsec)   : " << end - start << endl;
  cout << "MSE              : " << mse << endl;
  cout << "Repeat Count     : " << count << endl;
  cout << endl << "--------------" << endl;
  //end

  //for graph plot
  if(D == 2){
    ofstream ofs;
    ofs.open("output.txt");

    n = 0;
    while(n < N){
      ofs << x[n][0] << " " << x[n][1] << " " << S[n]+1 << endl;
      n++;
    }
  }
  //end
  
  return 0;
}

long double centroid(int degree, long double *trgt, long double *crtn){
  
  long double centroid = 0.0;
  
  int i = 0;
  while(i < degree){
    centroid += pow(trgt[i] - crtn[i],2);
    i++;
  }
  centroid = sqrt(centroid);
  
  return centroid;
}
