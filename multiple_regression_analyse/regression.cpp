/*

  Filename    :regression.cpp
  Filetype    :C++ Script
  Created     :20141125
  Updated     :20141125
  Author      :Reona SEKINO
  Copyright   :Copyright(C)2014 Reona SEKINO All Rights Reserved.
  Descliption :Regression Analysis Script for Large Scale Data
  Caution     :More improvement can be done in the imprementation.
  Format(IN)  :.csv(Comma-Separated Values)
  Format(OUT) :/txt(Plain-Text File)
  Usage       :$ ./a.out (Fill the name of Input CSV file)

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <math.h>
#include <sstream>

using namespace std;

int gauss(double **, int, int, double);
double *regression(int, int, double **, double *, double);

int main(int argc, char **argv)
{
  int i, j, h;
  
  /*** prepare for output file ***/
  ofstream ofs;
  ofs.open("result.txt");

  /*** read CSV file ***/
  vector<string> fieldnames;
  string buf;
  ifstream ifs;
  ifs.open(argv[1]);
  if(!ifs){
    cout << "Input data file not found" << endl;
    return 0;
  }
  
  /*** read fieldnames ***/
  {
    getline(ifs, buf);
    string token;
    istringstream stream(buf);
    
    while(getline(stream, token, ',')){
      fieldnames.push_back(token);
    }
  }

  /*** check field size ***/
  int n = fieldnames.size();
  vector< vector<double> > value;

  /*** read values ***/
  while(getline(ifs, buf)){
    double tmp;
    string token;
    istringstream stream(buf);

    vector<double> tmpArray;
    while(getline(stream, token, ',')){
      stringstream ss;
      ss << token;
      ss >> tmp;
      tmpArray.push_back(tmp);
    }
    value.push_back(tmpArray);
  }

  /*** check data size ***/
  int sizetmp = 0;
  int N = value.size();
  
  /*** installed data check ***/
  /*
    i = 0;
    while(i < N){
    if(n != value.at(i).size()){
    cout << "Format of input file is not supported" << endl;
    return 0;
    }
    i++;
    }
    i = 0;
    while(i < n){
    cout << fieldnames.at(i) << " ";
    i++;
    }
    cout << endl;
    i = 0;
    while(i < N){
    j = 0;
    while(j < n){
    cout << value.at(i).at(j) << " ";
    j++;
    }
    cout << endl;
    i++;
    }
  */

  /*** select response variable ***/
  int rspvar = -1;
  cout  << "Which data is used as Response Variable (equal to Y)? Please enter the field name to use as Response Variable." << endl;
  while(true){
    string checkstr;
    cin >> checkstr;
    i = 0;
    while(i < n){
      if(checkstr == fieldnames.at(i)){
	rspvar = i;
	cout << "Accepted." << endl;
	break;
      }
      i++;
    }
    if(rspvar == -1){
      cout << "There is no candidate for \"" << checkstr << "\". Please enter again." << endl;
    }
    else{
      break;
    }
  }

  /*** select explaination variable ***/
  vector<int> onfield;
  cout << "Please Answer if the field asked is to use or not." << endl;
  i = 0;
  while(i < n){
    char tmp;
    if(i == rspvar)i++;
    cout << "Do you use the \"" << fieldnames.at(i) << "\" field? (Please answer Y or N)" << endl;
    while(true){
      cin >> tmp;
      if(tmp == 'Y'){
	onfield.push_back(i);
	break;
      }
      else{
	if(tmp == 'N'){
	  break;
	}
	else{
	  cout << "Please answer \"Y\" or \"N\"." << endl;
	}
      }
    }
    i++;
  }

  /*** confirm the field ***/
  cout << "----------------" << endl;
  cout << "Response Variable : " << fieldnames.at(rspvar) << endl;
  cout << "Explaination Variable : ";
  i = 0;
  while(i < onfield.size()){
    cout << fieldnames.at(onfield.at(i)) << " ";
    i++;
  }
  cout << endl << "----------------" << endl;
  cout << "Analysis is started in above data field settings. Pleasr confirm in \"Y\" or \"N\"." << endl;
  while(true){
    char tmp;
    cin >> tmp;
    if(tmp == 'Y'){
      cout << "Now analysing......" << endl;
      break;
    }
    else{
      if(tmp == 'N'){
	cout << "Analysis is canceled. Please try again." << endl;
	return 0;
      }
      else{
	cout << "Please answer \"Y\" or \"N\"." << endl;
      }
    }
  }
  
  /*** calculate basical data ***/ 
  vector<double> stddev;
  vector<double> ave;
  
  /*** calculate average ***/
  i = 0;
  while(i < n){
    double avetmp = 0;
    j = 0;
    while(j < N){
      avetmp += value.at(j).at(i);
      j++;
    }
    ave.push_back(avetmp / N);
    i++;
  }

  /*** check average ***/
  /*
    i = 0;
    while(i < n){
    cout << "Average of " << fieldnames.at(i) << " : " << ave.at(i) << endl;
    i++;
    }
  */

  /*** calculate standard deviation ***/
  i = 0;
  while(i < n){
    double stddevtmp = 0;
    j = 0;
    while(j < N){
      stddevtmp += pow(value.at(j).at(i) - ave.at(i), 2);
      j++;
    }
    stddev.push_back(sqrt(stddevtmp / N));
    i++;
  }

  /*** check standard deviation ***/
  /*
    i = 0;
    while(i < n){
    cout << "Standard Deviation of " << fieldnames.at(i) << " : " << stddev.at(i) << endl;
    i++;
    }
  */

  ofs << "-----------------#RESULT#-----------------" << endl << endl;

  int count = 1;
  int flag = 1;
  while(flag){

    /*** regression with statdardization ***/
    {
      /*** initialize ***/
      n = onfield.size();
      double *b, *y, **X;
      y = new double [N];
      X = new double * [N];
      i = 0;
      while(i < N){
	X[i] = new double [n + 1];
	i++;
      }
    
      /*** set values ***/
      i = 0;
      while(i < N){
	X[i][0] = 1.0;
	y[i] = (value.at(i).at(rspvar) - ave.at(rspvar)) / stddev.at(rspvar); 
	j = 0;
	while(j < n){
	  X[i][j + 1] = (value.at(i).at(onfield.at(j)) - ave.at(onfield.at(j))) / stddev.at(onfield.at(j));
	  j++;
	}
	i++;
      }
    
      /*** check X and y***/
      /*
	i = 0;
	while(i < N){
	cout << y[i] << " ";
	j = 0;
	while(j < n){
	cout << X[i][j] << " ";
	j++;
	}
	cout << endl;
	i++;
	}
      */
    
      /*** calculate coefficient array ***/
      b = regression(n, N, X, y, 1.0e-10);   
    
      /*** check b ***/
      /*
	i = 0;
	while(i < n){
	cout << b[i] << endl;
	i++;
	}
      */

      if(b != NULL){//if regression is completed correctly
      
	/*** report information ***/
	ofs << "-------------<Main Informaiton>-----------" << endl << endl;
	ofs << "Count : " << count << endl;
	ofs << "Response Variable : " << fieldnames.at(rspvar) << endl;
	ofs << "Explaination Variable : ";
	i = 0;
	while(i < onfield.size()){
	  ofs << fieldnames.at(onfield.at(i)) << " ";
	  i++;
	}
	ofs << endl << "Number of field : " << n << endl;
	ofs << "Number of data set: " << N << endl;
	ofs << endl <<  "------------------------------------------" << endl << endl;
	ofs << "----------<With Standardization>----------" << endl << endl;
	ofs << "Constant Term = " << b[0] << endl;
	i = 0;
	while(i < n){
	  ofs << fieldnames.at(onfield.at(i)) << " = " << b[i + 1] << endl;
	  i++;
	}
	delete [] b;
      
      }
      else{//if regression is not completed (the data set is invalid)
     
	cout  << "Analysis is not completed because the loaded data set was not correct for multiple regression analysis." << endl;
	return 0;
      }

      /*** release memory ***/
      for (i = 0; i < N; i++)
	delete [] X[i];
      delete [] X;
      delete [] y;
    }

    ofs << endl <<  "------------------------------------------" << endl << endl;
  
    /*** regression without statdardization ***/
    {
      /*** initialize ***/
      n = onfield.size();
      double *b, *y, **X;
      y = new double [N];
      X = new double * [N];
      i = 0;
      while(i < N){
	X[i] = new double [n + 1];
	i++;
      }
    
      /*** set values ***/
      i = 0;
      while(i < N){
	X[i][0] = 1.0;
	y[i] = value.at(i).at(rspvar);
	j = 0;
	while(j < n){
	  X[i][j + 1] = value.at(i).at(onfield.at(j));
	  j++;
	}
	i++;
      }
    
      /*** check X and y***/
      /*
	i = 0;
	while(i < N){
	cout << y[i] << " ";
	j = 0;
	while(j < n){
	cout << X[i][j] << " ";
	j++;
	}
	cout << endl;
	i++;
	}
      */
    
      /*** calculate coefficient array ***/
      b = regression(n, N, X, y, 1.0e-10);   
    
      /*** check b ***/
      /*
	i = 0;
	while(i < n){
	cout << b[i] << endl;
	i++;
	}
      */

      if(b != NULL){//if regression is completed correctly

	/*** calculate t_value ***/
	double a[n + 1][n + 1];
	double Syy = 0;

	/*** initialize ***/
	i = 1;
	while(i < n + 1){
	  h = 1;
	  while(h < n + 1){
	    a[i][h] = 0;
	    j = 0;
	    while(j < N){
	      a[i][h] += (X[j][i] - ave.at(onfield.at(i - 1))) * (X[j][h] - ave.at(onfield.at(h - 1))); 	
	      j++;
	    }
	    h++;
	  }
	  i++;
	}

	/*** calculate reverse matrix ***/
	double beta[30];
	int ipv, j;
	double inv_pivot, temp;
	double big;
	int pivot_row, row[30];

	for(ipv=1 ; ipv <= n ; ipv++){

	  big=0.0;
	  for(i=ipv ; i<=n ; i++){
	    if(fabs(a[i][ipv]) > big){
	      big = fabs(a[i][ipv]);
	      pivot_row = i;
	    }
	  }
	  if(big == 0.0){return 0;}	
	  row[ipv] = pivot_row;

	  if(ipv != pivot_row){
	    for(i=1 ; i<=n ; i++){
	      temp = a[ipv][i];
	      a[ipv][i] = a[pivot_row][i];
	      a[pivot_row][i] = temp;
	    }
	    temp = beta[ipv];
	    beta[ipv] = beta[pivot_row];
	    beta[pivot_row] = temp;
	  }

	  inv_pivot = 1.0/a[ipv][ipv];
	  a[ipv][ipv]=1.0;
	  for(j=1 ; j <= n ; j++){
	    a[ipv][j] *= inv_pivot;
	  }
	  beta[ipv] *= inv_pivot;
	
	  for(i=1 ; i<=n ; i++){
	    if(i != ipv){
	      temp = a[i][ipv];
	      a[i][ipv]=0.0;
	      for(j=1 ; j<=n ; j++){
		a[i][j] -= temp*a[ipv][j];
	      }
	      beta[i] -= temp*beta[ipv];
	    }
	  }
	}
      
	for(j=n ; j>=1 ; j--){
	  if(j != row[j]){
	    for(i=1 ; i<=n ; i++){
	      temp = a[i][j];
	      a[i][j]=a[i][row[j]];
	      a[i][row[j]]=temp;
	    }
	  }
	}

	/*** calculate S_yy ***/
	i = 0;
	while(i < N){
	  double  tmp = 0.0;
	  j = 0;
	  while(j < n + 1){
	    tmp += X[i][j] * b[j];
	    j++;
	  }
	  Syy += pow(y[i] - tmp, 2);
	  i++;
	}
      
	/*** calculate coefficient of determination ***/
	double ssr = 0;
	double r2;
	double r2_fixed;
	i = 0;
	while(i < N){
	  double varprice = 0;
	  j = 0;
	  while(j < n + 1){
	    varprice += X[i][j] * b[j];
	    j++;
	  }
	  ssr += pow(y[i] - varprice, 2);
	  i++;
	}   
	r2 = 1 - (ssr / (pow(stddev[rspvar], 2) * N));
	r2_fixed = 1 - (ssr / (pow(stddev[rspvar], 2) * N)) * (N - 1) / (N - n - 1);

	/*** report information ***/
	ofs << "---------<Without standardization>--------" << endl << endl;
	ofs << "Constant Term = " << b[0] << endl;
	i = 0;
	while(i < n){
	  ofs << fieldnames.at(onfield.at(i)) << " = " << b[i + 1] << endl;
	  i++;
	}
	delete [] b;

	ofs << endl <<  "------------------------------------------" << endl << endl;
	ofs << "------<Coefficient of Determination>------" << endl << endl;
      

	ofs << "Degree of freedom unconsidered : " << r2 << endl;
	ofs << "Degree of freedom considered : " << r2_fixed << endl;
	ofs << endl <<  "------------------------------------------" << endl << endl;
	ofs << "--------<T Value and Significance>--------" << endl << endl;

	/*** calculate t_value ***/
	double t_value;
	vector<int> tmpfield;;
	flag = 0;
	i = 1;
	while(i < n + 1){
	  t_value = b[i] / sqrt(a[i][i] * Syy / (N - n - 1));
	  if(fabs(t_value) > 2.0){
	    ofs << fieldnames.at(onfield.at(i - 1)) << " : Significant（t_value = " << t_value << ")" << endl;
	    tmpfield.push_back(onfield.at(i - 1));
	  }
	  else{
	    ofs << fieldnames.at(onfield.at(i - 1)) << " : Unsignificant（t_value = " << t_value << ")" << endl;
	    flag = 1;
	  }
	  i++;
	}
	onfield.clear();
	onfield.swap(tmpfield);
	ofs << endl <<  "------------------------------------------" << endl << endl; 
      }
      else{//if regression is not completed (the data set is invalid)
     
	cout  << "Analysis is not completed because the loaded data set was not correct for multiple regression analysis." << endl;
	return 0;
      }

      /*** release memory ***/
      for (i = 0; i < N; i++)
	delete [] X[i];
      delete [] X;
      delete [] y;
    }
    ofs << endl;
    count ++;
  }
  cout << "Analysis is completed. Please check the output file (result.txt)." << endl;
  ofs << "Analysed by regression.cpp" << endl;
  ofs << "Copyright(C) 2014 Reona SEKINO All Right Reserved." << endl;

  return 0;
}


double *regression(int n, int N, double **X, double *y, double eps)
{
  double **w, *b;
  int i1, i2, i3, sw;

  n++;
  b = new double [n];
  w = new double * [n];
  for (i1 = 0; i1 < n; i1++)
    w[i1] = new double [n+1];

  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 0; i2 < n; i2++) {
      w[i1][i2] = 0.0;
      for (i3 = 0; i3 < N; i3++)
	w[i1][i2] += X[i3][i1] * X[i3][i2];
    }
  }

  for (i1 = 0; i1 < n; i1++) {
    w[i1][n] = 0.0;
    for (i2 = 0; i2 < N; i2++)
      w[i1][n] += X[i2][i1] * y[i2];
  }
  
  sw = gauss(w, n, 1, eps);

  if (sw == 0) {
    for (i1 = 0; i1 < n; i1++)
      b[i1] = w[i1][n];
  }
  else
    b = NULL;

  for (i1 = 0; i1 < n; i1++)
    delete [] w[i1];
  delete [] w;

  return b;
}
  
int gauss(double **w, int n, int m, double eps)
{
  double y1, y2;
  int ind = 0, nm, m1, m2, i1, i2, i3;

  nm = n + m;

  for (i1 = 0; i1 < n && ind == 0; i1++) {

    y1 = .0;
    m1 = i1 + 1;
    m2 = 0;

    for (i2 = i1; i2 < n; i2++) {
      y2 = fabs(w[i2][i1]);
      if (y1 < y2) {
	y1 = y2;
	m2 = i2;
      }
    }

    if (y1 < eps)
      ind = 1;

    else {

      for (i2 = i1; i2 < nm; i2++) {
	y1        = w[i1][i2];
	w[i1][i2] = w[m2][i2];
	w[m2][i2] = y1;
      }

      y1 = 1.0 / w[i1][i1];

      for (i2 = m1; i2 < nm; i2++)
	w[i1][i2] *= y1;

      for (i2 = 0; i2 < n; i2++) {
	if (i2 != i1) {
	  for (i3 = m1; i3 < nm; i3++)
	    w[i2][i3] -= w[i2][i1] * w[i1][i3];
	}
      }
    }
  }

  return(ind);
}
