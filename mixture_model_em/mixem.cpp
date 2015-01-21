#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

using namespace std;

int main(void)
{
  ifstream ifs;
  ifs.open("./sample_mix.txt");

  ofstream ofs;
  ofs.open("./output.txt");
  
  char trash;
  ifs >> trash;
  ifs >> trash;

  vector<long double> x;
  vector<long double> y;

  while(!(ifs.eof())){
    long double temp;
    ifs >> temp;
    x.push_back(temp);
    ifs >> temp;
    y.push_back(temp);
  }

  int N = x.size() - 1;//data set size(200)
  int K = 5;//size of cross validation
  int E = N / K;//for cross validation(20)

  long double L[N / 10];//log likelihood
  long double H[N / 10][E];//Cross Entorpy(Ave.)
  long double paramArray[N / 10][E][N / 10][6];//save parameters

  ifs.close();//installing sample data
  
  int M = 1;
  while(M < N / 10 + 1){
    int e = 0;
    while(e < K){

      long double mu_x[M];
      long double mu_y[M];//average
      long double sigma_xx[M];
      long double sigma_xy[M];
      long double sigma_yy[M];//variance-covariance matrix factor
      long double pi[M];
      long double pi_temp[M];//temporary
      long double delta;//determinant
      long double f[N][M];//Fn
      long double w[M][N];

      int m = 1;
      int n = 1;
      long double mc = 1;
      while(m < M + 1){
	mu_x[m - 1] = (2 * mc - 1) / (2 * M);
	mu_y[m - 1] = (2 * mc - 1) / (2 * M);
	m++; 
	mc++;
      }//format mu

      m = 1;
      while(m < M + 1){
	sigma_xx[m - 1] = 0;
	sigma_xy[m - 1] = 0;
	sigma_yy[m - 1] = 0;
	n = 1;
	while(n < N + 1){
	  if(n % K != e){
	    sigma_xx[m - 1] += (x[n - 1] - mu_x[m - 1]) * (x[n - 1] - mu_x[m - 1]);
	    sigma_xy[m - 1] += (x[n - 1] - mu_x[m - 1]) * (y[n - 1] - mu_y[m - 1]);
	    sigma_yy[m - 1] += (y[n - 1] - mu_y[m - 1]) * (y[n - 1] - mu_y[m - 1]);
	  }
	  n++;
	}
	sigma_xx[m - 1] = sigma_xx[m - 1] / (N - E);
	sigma_xy[m - 1] = sigma_xy[m - 1] / (N - E);
	sigma_yy[m - 1] = sigma_yy[m - 1] / (N - E);
	
	m++;
      }//format sigma

      m = 1;
      while(m < M + 1){
	pi[m - 1] = 1.0 / M;
	m++;
      }//format pi

      long double L_temp = 0;
      int fflag = 1;
      int count = 0;
      m = 1;
      while(m < M + 1){
	if((M == 3) && (e == 4))
	  ofs << mu_x[m - 1] << " " << mu_y[m - 1] << " " << sigma_xx[m - 1] << " " << sigma_yy[m - 1] << " " << sigma_xy[m - 1] << " " << pi[m - 1] << endl;
	m++;
      }

      while(1){//convergence requirement
	count++;
	
	m = 1;
	while(m < M + 1){
	  delta = sigma_xx[m - 1] * sigma_yy[m - 1] - sigma_xy[m - 1] * sigma_xy[m - 1];
	  n = 1;
	  while(n < N + 1){
	    if(n % K != e){
	      f[n - 1][m - 1] = exp((2 * sigma_xy[m - 1] * (x[n - 1] - mu_x[m - 1]) * (y[n - 1] - mu_y[m - 1]) - sigma_xx[m - 1] * (y[n - 1] - mu_y[m - 1]) * (y[n - 1] - mu_y[m - 1]) - sigma_yy[m - 1] * (x[n - 1] - mu_x[m - 1]) * (x[n - 1] - mu_x[m - 1])) / (2 * delta)) / (2 * M_PI * sqrt(delta));
	    }
	    n++;
	  }
	  m++;
	}//calculate Fn
     
	n = 1;
	while(n < N + 1){
	  if(n % K != e){
	    long double sgm = 0;
	    m = 1;
	    while(m < M + 1){
	      sgm += f[n - 1][m - 1] * pi[m - 1];
	      m++;
	    }
	    m = 1;
	    while(m < M + 1){
	      w[m - 1][n - 1] = (f[n - 1][m - 1] * pi[m - 1]) / sgm;
	      m++;
	    }
	  }
	  n++;
	}//calculate w //E step

	m = 1;
	while(m < M + 1){
	  long double den = 0;
	  n = 1;
	  while(n < N + 1){
	    if(n % K != e){
	      den += w[m - 1][n - 1];
	    }
	    n++;
	  }//calculate denominator

	  mu_x[m - 1] = 0;
	  mu_y[m - 1] = 0;

	  n = 1;
	  while(n < N + 1){
	    if(n % K != e){
	      mu_x[m - 1] += (x[n - 1] * w[m - 1][n - 1]) / den;
	      mu_y[m - 1] += (y[n - 1] * w[m - 1][n - 1]) / den;//update mu
	    }
	    n++;
	  }

	  sigma_xx[m - 1] = 0;
	  sigma_xy[m - 1] = 0;
	  sigma_yy[m - 1] = 0;

	  n = 1;
	  while(n < N + 1){
	    if(n % K != e){
	      sigma_xx[m - 1] += ((x[n - 1] - mu_x[m - 1]) * (x[n - 1] - mu_x[m - 1]) * w[m - 1][n - 1]) / den;
	      sigma_xy[m - 1] += ((x[n - 1] - mu_x[m - 1]) * (y[n - 1] - mu_y[m - 1]) * w[m - 1][n - 1]) / den;
	      sigma_yy[m - 1] += ((y[n - 1] - mu_y[m - 1]) * (y[n - 1] - mu_y[m - 1]) * w[m - 1][n - 1]) / den;
	    }
	    n++;
	  }//update sigma

	  pi[m - 1] = den / (N - E);//update pi
	  m++;
	}//M step

	m = 1;
	while(m < M + 1){
	  if((M == 3) && (e == 4))
	    ofs << mu_x[m - 1] << " " << mu_y[m - 1] << " " << sigma_xx[m - 1] << " " << sigma_yy[m - 1] << " " << sigma_xy[m - 1] << " " << pi[m - 1] << endl;
	  m++;
	}
     
	L[M - 1] = 0;
	n = 1;
	while(n < N + 1){
	  if(n % K == e){
	    long double sgm = 0;
	    m = 1;
	    while(m < M + 1){
	      delta = sigma_xx[m - 1] * sigma_yy[m - 1] - sigma_xy[m - 1] * sigma_xy[m - 1];
	      sgm += pi[m - 1] * exp((2 * sigma_xy[m - 1] * (x[n - 1] - mu_x[m - 1]) * (y[n - 1] - mu_y[m - 1]) - sigma_xx[m - 1] * (y[n - 1] - mu_y[m - 1]) * (y[n - 1] - mu_y[m - 1]) - sigma_yy[m - 1] * (x[n - 1] - mu_x[m - 1]) * (x[n - 1] - mu_x[m - 1])) / (2 * delta)) / (2 * M_PI * sqrt(delta));
	      m++;
	    }
	    L[M - 1] += log(sgm);
	  }
	  n++;
	}//calculate L

	if(fflag){
	  L_temp = L[M - 1];
	  fflag = 0;
	} 
	else{
	  if(L[M - 1] - L_temp < 0){
	    L[M - 1] = L_temp;
	    break;
	  }
	  if(L[M - 1] - L_temp > 0.01){//termination judgement
	    L_temp = L[M - 1];
	  }
	  else{
	    break;
	  }
	}
      
      }//complete maximize L

      m = 1;
      while(m < M + 1){
	paramArray[M - 1][e][m - 1][0] = mu_x[m - 1];
	paramArray[M - 1][e][m - 1][1] = mu_y[m - 1];
	paramArray[M - 1][e][m - 1][2] = sigma_xx[m - 1];
	paramArray[M - 1][e][m - 1][3] = sigma_yy[m - 1];
	paramArray[M - 1][e][m - 1][4] = sigma_xy[m - 1];
	paramArray[M - 1][e][m - 1][5] = pi[m - 1];
	m++;
      }
      H[M - 1][e] =  - L[M - 1] / E;//caluculate cross entropy
      //cout << M << " " << e << " " << H[M - 1][e] << endl;
      e++;
    }//complete cross validation

    M++;
  }

  int i = 1;
  long double min = 0xffffffff;
  M = 1;
  while(M < N / 10 + 1){
    long double tmp = 0.0;
    int e = 1;
    while(e < E){
      tmp += H[M - 1][e];
      e++;
    }
    if(tmp / E < min){
      min = tmp / E;
      i = M;
    }
    M++;
  }//most fitting quantity of factor

  int j = 0;
  int e = 1;
  min = H[i - 1][0];
  while(e < E){
    if(H[i - 1][e] < min){
      j = e;
      min = H[i - 1][e];
    }
    e++;
  }

  cout << "quantity of factor           : " << i << "         (used e = " << j << ")" << endl;
  int k = 0;
  while(k < i){
    cout << "Factor " << k + 1 << endl;
    cout << "(mu_x,mu_y)                  : (" << paramArray[i - 1][j][k][0] << "," << paramArray[i - 1][j][k][1] << ")" << endl;
    cout << "(sigma_xx,sigma_yy,sigma_xy) : (" << paramArray[i - 1][j][k][2] << "," << paramArray[i - 1][j][k][3] << "," << paramArray[i - 1][j][k][4] << ")" << endl;
    cout << "pi                           : " << paramArray[i - 1][j][k][5] << endl;
    k++;
  }
  cout<<  "likelihood function          : " << L[i - 1] << " " << endl;

  return 0;
}
