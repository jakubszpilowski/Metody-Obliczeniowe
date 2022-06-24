#include "func.h"

double **get_matrix_conditions(int n, int m){
	double **matrix = new double *[n];
	
	for(int i = 0; i < n; i++){
		matrix[i] = new double[m];
	}
		
	/* Warunki początkowe do połowy kolumn w macierzy 1.0, potem 0.0
	* U(x,0) = 1.0 dla x < 0
	* U(x,0) = 0.0 dla x >= 0
	*/
	int middle = m / 2;
	for(int i = 0; i < middle; i++){
		matrix[0][i] = 1.0;
	}
	
	for(int i = middle; i < m; i++){
		matrix[0][i] = 0.0;
	}
	

	/* Warunki brzegowe, na końcach przedziałów 1.0 i 0.0
	* U(-a, t) = 1.0
	* U(a, t) = 0.0
	*/
	for(int i = 0; i < n; i++){
		matrix[i][0] = 1.0;
		matrix[i][m - 1] = 0.0;
	}
	
	return matrix;
}


double get_delta_t(const double lambda, const double h, const double D){
	//Przekształcamy wzór LAMBDA = D * delta_t / (h * h)
	return (lambda * h * h) / D;
}


double **get_matrix(int n, int m){
	double **matrix = new double *[n];
	
	for(int i = 0; i < n; i++){
		matrix[i] = new double[m];
	}
	
	return matrix;
}


double *get_vector(int m){
	double *x = new double [m];
	
	for(int i = 0; i < m; i++)
		x[i] = 0.0;
	
	return x;
}


double **get_error_matrix(double **analytical, double **result, const int n, const int m){
	double **err = get_matrix(n, m);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			err[i][j] = fabs(analytical[i][j] - result[i][j]);
		}
	}
	
	return err;
}


double *get_error_vector(double **err, const int n, const int m){
	double *v = new double[n];
	
	for(int i = 0; i < n; i++){
		double current_max = fabs(err[i][0]);
		for(int j = 1; j < m; j++){
			current_max = max(current_max, fabs(err[i][j]));
		}
		v[i] = current_max;
	}
	
	return v;
}


double *error_vector(double *a, double *x, int m){
	double *v = new double[m];

	for(int i = 0; i < m; i++)
		v[i] = fabs(a[i] - x[i]);

	return v;
}


double get_max_error(double *v, const int m){
	double max = fabs(v[0]);

	for(int i = 1; i < m; i++){
		if(fabs(v[i]) > max){
			max = fabs(v[i]);
		}
	}

	return log10(max);
}


double *get_steps_t(double delta_t, int n, double t_0){
	double *steps = new double[n];
	
	double t = t_0;
	
	for(int i = 0; i < n; i++){
		steps[i] = t;
		t += delta_t;
	}
	
	return steps;
}


void clean_up(double **results, double **a, double **b, double **c, double **d, int n){
	for(int i = 0; i < n; i++){
		delete [] results[i];
		delete [] a[i];
		delete [] b[i];
		delete [] c[i];
		delete [] d[i];
	}
	delete [] results;
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] d;
}


void delete_matrix(double **matrix, int n){
	for(int i = 0; i < n; i++)
		delete[] matrix[i];
	delete[] matrix;
}


void write_results(double *analitical, double *numerical, double *err, int m, double h, double x_start, const char* filename){
	fstream output(filename, ios::out);
	
	double x = x_start;
	for(int i = 0; i < m ; i++){
		output << setprecision(10) << x << " " << analitical[i];
		output << " " << numerical[i] << " " << err[i] << endl;
		x += h;
	}
	output.close();
}


bool check_en(double *x_n, double* x_n1, int m){
	int counter = 0;
	
	for(int i = 0; i < m; i++){
		if(fabs(x_n1[i] - x_n[i]) < TOLX)
			counter++;
	}

	return counter == m;
}


bool check_residue(double* l, double *d, double *u, double *x, double* b, int m){
	int counter = 0;
	
	//Zamiast używać całej macierzy m x m korzytsam z 3 wektorów
	//Pierwsze mnożenie nie zawiera elementu z wektora l dlatego pomijamy
	double tmp = d[0] * x[0] + u[0] * x[1];
	if(fabs(tmp - b[0]) < TOLF)
			counter++;

	//Mnożenie wektora rozwiązania, za każdym razem mnożymy tylko 3 wartości
	for(int i = 1; i < m - 1; i++){
		tmp = l[i] * x[i - 1] + d[i] * x[i] + u[i] * x[i + 1];
		if(fabs(tmp - b[i]) < TOLF)
			counter++;
	}
	
	//Ostatnie mnożenie nie zawiera elementu wektora u, pomijamy
	tmp = d[m - 1] * x[m - 1] + l[m - 1] * x[m - 2];
	if(fabs(tmp - b[m - 1]) < TOLF)
			counter++;
	
	return counter == m;
}