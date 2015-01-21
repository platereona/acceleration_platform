#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

using namespace std;

int main(void)
{
  int n = 0;
  int m = 0;
  int M = 3;
  int N = 200;//data set size(200)
  long double mu_x[M];
  long double mu_y[M];//average
  long double sigma_xx[M];
  long double sigma_xy[M];
  long double sigma_yy[M];//variance-covariance matrix factor
  long double pi[M];
  long double P[N][M];
  int count[M];
  
  
  ifstream ifs;
  ifs.open("./sample_mix.txt");

  ifstream ifs2;
  ifs2.open("./param.txt");

  ofstream ofs;
  ofs.open("./result.txt");
  
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

  m = 1;
  while(m < M + 1){
    ifs2 >> mu_x[m - 1];
    ifs2 >> mu_y[m - 1];
    ifs2 >> sigma_xx[m - 1];
    ifs2 >> sigma_yy[m - 1];
    ifs2 >> sigma_xy[m - 1];
    ifs2 >> pi[m - 1];
    m++;
  }

  n = 1;
  while(n < N + 1){
    m = 1;
    while(m < M + 1){
      long double delta;
      delta = sigma_xx[m - 1] * sigma_yy[m - 1] - sigma_xy[m - 1] * sigma_xy[m - 1];
      P[n - 1][m - 1] = pi[m - 1] * exp((2 * sigma_xy[m - 1] * (x[n - 1] - mu_x[m - 1]) * (y[n - 1] - mu_y[m - 1]) - sigma_xx[m - 1] * (y[n - 1] - mu_y[m - 1]) * (y[n - 1] - mu_y[m - 1]) - sigma_yy[m - 1] * (x[n - 1] - mu_x[m - 1]) * (x[n - 1] - mu_x[m - 1])) / (2 * delta)) / (2 * M_PI * sqrt(delta));
      m++;
    }

    long double max = 0.0;
    int maxtag = 0;
    m = 1;
    while(m < M + 1){
      if(P[n - 1][m - 1] > max){
	max = P[n - 1][m - 1];
	maxtag = m;
      }
      m++;
    }

    count[maxtag - 1]++;
    ofs << x[n - 1] << " " << y[n - 1] << " " << maxtag << endl;
  
    n++;
  }

  m = 1;
  while(m < M + 1){
    cout << m << " : " << count[m - 1] << endl;
    m++;
  }

  return 0;
}
