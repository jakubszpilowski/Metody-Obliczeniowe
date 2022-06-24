#include <iostream>
#include <cmath>
#include <fstream>
#include "calerf.h" 
#include "func.h"

using namespace std;

//Stałe
#define MAX_ITER 50 //maksymalna ilość iteracji dla gaussa - seidla

//współczynnik transportu z treści zadania
const double D = 1.0;

//przedział czasowy z treści zadania
const double t_0 = 0.0;
const double t_max = 2.0;

//przedział przestrzenny z treści zadania
const double a = 6.0 * sqrt(D * t_max);
const double x_start = -a;
const double x_end = a;

//współczynnik lambda dla metod pośrednich
const double LAMBDA = 1.0;


/**
 * Funkcja obliczająca macierz wyników analaitycznych
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn macierzy
 * @param h krok przestrzenny
 * @param delta_t krok czasowy
 * @return double** macierz wynikowa
 */
double **get_analytical_solution(const int n, const int m, const double h, const double delta_t){
	double **matrix = get_matrix(n, m);
	
	double x = x_start;
	double t = t_0;
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			//Obliczamy wartość korzystając z funkcji ERFC_L z pakietu CALERF 
			matrix[i][j] = 0.5 * calerf::ERFC_L(x / (2.0 * sqrt(D * t)));
			x += h;
		}
		x = x_start;
		t += delta_t;
	}
	
	return matrix;
}


/*=================================
*	ALGORYTMY
* =================================
*/

/**
 * Algorytm Thomasa z laboratorium numer 6, korzystam z 3 wektorów o rozmiarze m
 * @param lower wartości pod przekątną główną
 * @param diagonal wartości na przekątnej głownej
 * @param upper wartości nad przekątną głowną
 * @param b wetor wyniku
 * @param x wektor rozwiązania
 * @param m rozmiar wektora
 */
void Thomas(double *lower, double *diagonal, double *upper, double *b, double *x, int m){
	//Wyliczamy wartości eta na przekątnej głównej
	for(int i = 1; i < m; i++){
		diagonal[i] = diagonal[i] - (upper[i - 1] * lower[i]) / diagonal[i - 1];
	}
	
	//Wyliczamy współczynniki r w wektorze b
	for(int i = 1; i < m; i++){
		b[i] = b[i] - (lower[i] * b[i - 1]) / diagonal[i - 1];
	}
	
	//Wyliczamy x
	x[m - 1] = b[m - 1] / diagonal[m - 1];
	for(int i = m - 2; i >= 0; i--){
		x[i] = (b[i] - upper[i] * x[i + 1]) / diagonal[i];
	}
}


/**
 *	Algorytm iteracyjny Gaussa - Seidela, korzystam z 3 wektorów zamiast całej macierzy
 *	wersja z (L + D)x + Ux (analogicznie można na Lx + (D + U)x)
 * @param lower wartości pod przekątną główną
 * @param diagonal wartości na przekątnej głównej
 * @param upper wartości nad przekątną główną
 * @param b wektor wyniku
 * @param x wektor rozwiązania
 * @param m rozmiar wektorów
 */
void Gauss_Seidel_method(double *lower, double *diagonal, double *upper, double *b, double *x, int m){
	//Wektor wartości po prawej stronie równania
	double *right_side = get_vector(m);
	//Wektor kolejnego przybliżenia
	double *x_n1 = get_vector(m);
	
	for(int i = 0; i < MAX_ITER; i++){
		//Obliczamy wartości -U * x_n
		for(int j = 0; j < m - 1; j++){
			right_side[j] = -upper[j] * x[j + 1];
		}
		
		//Ostatni wyraz z prawej strony nie jest mnożony przez x
		right_side[m - 1] = -upper[m - 1];
		
		//Obliczmy wartości po prawej stronie równanie
		for(int j = 0; j < m; j++)
			right_side[j] += b[j];
			
		//Rozwiązujemy układ z macierzą trójkątną dolną
		x_n1[0] = right_side[0] / diagonal[0];
			
		for(int j = 1; j < m; j++){
			x_n1[j] = (right_side[j] - lower[j] * x_n1[j - 1]) / diagonal[j];
		}
		
		//Sprawdzamy kryterium dokładności
		bool e_n = check_en(x, x_n1, m);

		//Sprawdzamy kryterium wiarygodności
		bool residue = check_residue(lower, diagonal, upper, x, b, m);
		
		//Jeśli w normie wychodzimy z metody - znaleźliśmy dobre przybliżenie
		if(e_n && residue){
			break;
		}

		//Jeśli nie wektor x_n1 zamieniamy z x i na nim wykonujemy obliczenia	
		for(int i = 0; i < m; i++)
			swap(x[i], x_n1[i]);
	}

	//Usuwamy zaalokowane wektory pomocnicze
	delete [] right_side;
	delete [] x_n1;
}


/* =========================
*	DYSKRETYZACJA
*  =========================
*/
   
/**
 * Funkcja obliczająca numerycznie wartości równania różniczkowego 
 * korzystając z algorytmu Thomasa i dyskretyzacji metodą Laasonen
 * @param n liczba wierszy macierzy (ilość kroków czasowych)
 * @param m liczba kolumn macierzy (ilość kroków przestrzennych)
 * @return double** macierz wyników numerycznych
 */
double **get_Laasonen_Thomas(int n, int m){
	//alokujemy macierz z warunkami brzegowymi i początkowymi
	double **matrix_U = get_matrix_conditions(n, m);
	
	//w celu uniknięcia liczenia tego wyrażenia w każdej iteracji liczymy tylko raz przed pętlą
	const double lambda = 1.0 + 2.0 * LAMBDA;
	
	//alokujemy potrzebne wektory
	double *lower = new double[m];
	double *diagonal = new double[m];
	double *upper = new double[m];
	double *b = new double[m];
	double *x = new double[m];
	
	//idziemy od poziomu t = t_0 + delta_t <- t = t_0 już mamy
	for(int k = 1; k < n; k++){
		//wypełniamy poszczególne wektory zgodnie z warunkami brzegowymi
		lower[0] = 0.0;
		diagonal[0] = 1.0;
		upper[0] = 0.0;
		b[0] = matrix_U[k - 1][0];
		
		for(int i = 1; i < m - 1; i++){
			//wypełniamy wektory zgodnie z zasadą dyskretyzacji metody Laasonen
			lower[i] = LAMBDA;
			diagonal[i] = -lambda;
			upper[i] = LAMBDA;
			b[i] = -matrix_U[k - 1][i];
		}
		
		//Wypełniamy zgodnie z warunkami brzegowymi
		lower[m - 1] = 0.0;
		diagonal[m - 1] = 1.0;
		upper[m - 1] = 0.0;
		b[m - 1] = matrix_U[k - 1][m - 1];
		
		//obliczamy rozwiązanie dla poziomu czasowego t = t_0 + k * delta_t
		Thomas(lower, diagonal, upper, b, x, m);
		
		//obliczone wartości zapisujemy w odpowiednie miejsce w macierzy wynikowej
		for(int i = 1; i < m - 1; i++)
			matrix_U[k][i] = x[i];
	}
	
	//usuwamy zaalokowane wektory
	delete [] lower;
	delete [] upper;
	delete [] diagonal;
	delete [] b;
	delete [] x;

	//zwracamy macierz z obliczonymi wynikami
	return matrix_U;
}


/**
 * Funkcja obliczająca numerycznie wartości równania różniczkowego 
 * korzystając z algorytmu Thomasa i dyskretyzacji metodą Cranka - Nicolson
 * @param n liczba kroków czasowych
 * @param m liczba kroków przestrzennych
 * @return double** macierz wyników numerycznych
 */
double **get_Crank_Nicolson_Thomas(int n, int m){
	//alokujemy macierz z odpowiednimi warunkami
	double **matrix_U = get_matrix_conditions(n, m);
	
	//obliczymy raz wartości współczynników pomocniczych
	const double lambda_half = LAMBDA / 2.0;
	const double lambda_plus = 1.0 + LAMBDA;
 	
	//alokujemy potrzebne wektory
	double *lower = new double[m];
	double *diagonal = new double[m];
	double *upper = new double[m];
	double *b = new double[m];
	double *x = new double[m];
	
	//podobnie jak dla Laasonen zaczynamy od poziomu czasowego t_0 + delta_t
	for(int k = 1; k < n; k++){
		//wypełniamy wektory odpowiednio z warunku brzegowego
		lower[0] = 0.0;
		diagonal[0] = 1.0;
		upper[0] = 0.0;
		b[0] = matrix_U[k - 1][0];
		
		for(int i = 1; i < m - 1; i++){
			//wypełniamy zgodnie z zasadami dyskretyzacji metody Cranka-Nicolson
			lower[i] = lambda_half;
			diagonal[i] = -lambda_plus;
			upper[i] = lambda_half;
			//współczynnik (1 - LAMBDA) się zeruje dlatego możemy usunąć U[k - 1][i]
			b[i] = - (lambda_half * matrix_U[k - 1][i - 1] + lambda_half * matrix_U[k - 1][i + 1]);
		}
		
		//wypełniamy z drugiego warunku brzegowego
		lower[m - 1] = 0.0;
		diagonal[m - 1] = 1.0;
		upper[m - 1] = 0.0;
		b[m - 1] = matrix_U[k - 1][m - 1];
		
		//obliczamy rozwiązanie dla danego poziomu czasowego
		Thomas(lower, diagonal, upper, b, x, m);
		
		//obliczone wartości zapisujemy w odpowiednim miejscu macierzy
		for(int i = 1; i < m - 1; i++)
			matrix_U[k][i] = x[i];
	}
	
	//usuwamy wektory
	delete [] lower;
	delete [] diagonal;
	delete [] upper;
	delete [] b;
	delete [] x;

	//zwracamy macierz wyników
	return matrix_U;
}


/**
 * Funkcja obliczająca numerycznie wartości równania różniczkowego 
 * korzystając z algorytmu iteracyjnego Gaussa - Seidela i dyskretyzacji metodą Laasonen
 * @param n ilość kroków czasowych
 * @param m ilość kroków przestrzennych
 * @return double** macierz wyników numerycznych
 */
double** get_Laasonen_Gauss_Seidel(int n, int m){
	//alokujemy macierz początkową z warunkami
	double **matrix_U = get_matrix_conditions(n, m);
	
	//obliczamy pomocniczy współczynnik
	const double lambda = 1.0 + 2.0 * LAMBDA;
	
	//alokujemy potrzebne wektory
	double *lower = new double[m];
	double *diagonal = new double[m];
	double *upper = new double[m];
	double *b = new double[m];
	double *x = new double[m];
	
	//wektor x inicjujemy wartościami z kroku czasowego t_0 jako przybliżenie początkowe
	for(int i = 0; i < m; i++)
		x[i] = matrix_U[0][i];
	
	//obliczamy wartości zaczynając z kroku czasowego t_0 + delta_t
	for(int k = 1; k < n; k++){
		//wypełniamy wektory według warunku brzegowego
		lower[0] = 0.0;
		diagonal[0] = 1.0;
		upper[0] = 0.0;
		b[0] = matrix_U[k - 1][0];
		
		for(int i = 1; i < m - 1; i++){
			//wypełniamy zgodnie z dyskretyzacją metody Laasonen
			lower[i] = LAMBDA;
			diagonal[i] = -lambda;
			upper[i] = LAMBDA;
			b[i] = -matrix_U[k - 1][i];
		}
		
		//wypełniamy zgodnie z warunkiem brzegowym
		lower[m - 1] = 0.0;
		diagonal[m - 1] = 1.0;
		upper[m - 1] = 0.0;
		b[m - 1] = matrix_U[k - 1][m - 1];
		
		//obliczamy wartości równania
		Gauss_Seidel_method(lower, diagonal, upper, b, x, m);
		
		//wyniki zapisujemy w macierzy
		for(int i = 1; i < m - 1; i++)
			matrix_U[k][i] = x[i];
	}

	//usuwamy wektory
	delete [] lower;
	delete [] diagonal;
	delete [] upper;
	delete [] b;
	delete [] x;
	
	//zwracamy macierz wynikową
	return matrix_U;
}


/**
 * Funkcja obliczająca numerycznie wartości równania różniczkowego 
 * korzystając z algorytmu iteracyjnego Gaussa - Seidela i dyskretyzacji metodą Cranka - Nicolson
 * @param n ilość kroków czasowych
 * @param m ilość kroków przestrzennych
 * @return double** macierz wyników numerycznych
 */
double **get_Crank_Nicolson_Gauss_Seidel(int n, int m){
	//pobieramy macierz początkową
	double **matrix_U = get_matrix_conditions(n, m);
	
	//obliczamy współczynniki pomocnicze
	const double lambda_half = LAMBDA / 2.0;
	const double lambda_plus = 1.0 + LAMBDA;
 	
	//alokujemy potrzebne wektory
	double *lower = new double[m];
	double *diagonal = new double[m];
	double *upper = new double[m];
	double *b = new double[m];
	double *x = new double[m];
	
	//pobieramy z macierzy przybliżenie początkowe
	for(int i = 0; i < m; i++)
		x[i] = matrix_U[0][i];
	
	//obliczamy wartości zaczynając od poziomu t_0 + delta_t
	for(int k = 1; k < n; k++){
		//wypełniamy wektory według warunku brzegowego
		lower[0] = 0.0;
		diagonal[0] = 1.0;
		upper[0] = 0.0;
		b[0] = matrix_U[k - 1][0];
		
		for(int i = 1; i < m - 1; i++){
			//wypełniamy zgodnie z zasadami dyskretyzacji metody Cranka-Nicolson
			lower[i] = lambda_half;
			diagonal[i] = -lambda_plus;
			upper[i] = lambda_half;
			//podobnie jak dla Thomasa (1 - LAMBDA) się zeruje - usuwamy U[k - 1][i]
			b[i] = - (lambda_half * matrix_U[k - 1][i - 1] + lambda_half * matrix_U[k - 1][i + 1]);
		}
		
		//wypełniamy z drugiego warunku brzegowego
		lower[m - 1] = 0.0;
		diagonal[m - 1] = 1.0;
		upper[m - 1] = 0.0;
		b[m - 1] = matrix_U[k - 1][m - 1];
		
		//obliczamy wartości dla danego poziomu czasowego
		Gauss_Seidel_method(lower, diagonal, upper, b, x, m);
		
		//wyniki zapisujemy w macierzy
		for(int i = 1; i < m - 1; i++)
			matrix_U[k][i] = x[i];
	}

	//usuwamy wektory
	delete [] lower;
	delete [] diagonal;
	delete [] upper;
	delete [] b;
	delete [] x;
	
	//zwracamy macierz wynikową
	return matrix_U;
}


/* =====================
* 	OBLICZANIE I ZAPIS DO PLIKÓW
*  =====================
*/


/**
 * Funkcja zapisująca do pliku wartości analityczne i numeryczne dla różnych 
 * chwili t z całego przedziału
 * @param analytical macierzy wyników analitycznych
 * @param thomas_l macierzy wyników numerycznych dla alg Thomasa i dyskretyzacji Laasonen
 * @param thomas_cn macierz wyników numerycznych dla alg Thomasa i dyskretyzacji Cranka-Nicolson
 * @param gauss_l macierzy wyników numerycznych dla metody Gaussa-Seidela i dyskretyzacji Laasonen
 * @param gauss_cn macierzy wyników numerycznych dla metody Gaussa-Seidela i dyskretyzacji Cranka-Nicolson
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn macierzy
 * @param h krok h
 */
void draw_plot_time(double **analytical, double **thomas_l, double **thomas_cn, double **gauss_l, double **gauss_cn, int n, int m, double h){
	fstream output;
	//Do pliku zapisujemy wartości z różnych poziomów czasowych

	// t = 0.2s
	double x = x_start;
	int t = 0.1 * (n - 1);
	output.open("wykres_t02.txt", ios::out);
	for(int i = 0; i < m; i+=5){
		output << setprecision(8) << x << " " << analytical[t][i] << " ";
		output << thomas_l[t][i] << " " << thomas_cn[t][i] << " ";
		output << gauss_l[t][i] << " " << gauss_cn[t][i] << endl;
		x += h;
	}
	output.close();

	// t = 0.4s
	x = x_start;
	t = 0.2 * (n - 1);
	output.open("wykres_t04.txt", ios::out);
	for(int i = 0; i < m; i+=5){
		output << x << " " << analytical[t][i] << " ";
		output << thomas_l[t][i] << " " << thomas_cn[t][i] << " ";
		output << gauss_l[t][i] << " " << gauss_cn[t][i] << endl;
		x += h;
	}
	output.close();

	// t = 1s
	x = x_start;
	t = 0.5 * (n - 1);
	output.open("wykres_t1.txt", ios::out);
	for(int i = 0; i < m; i+=5){
		output << setprecision(8) << x << " " << analytical[t][i] << " ";
		output << thomas_l[t][i] << " " << thomas_cn[t][i] << " ";
		output << gauss_l[t][i] << " " << gauss_cn[t][i] << endl;
		x += h;
	}
	output.close();

	// t = 1.6s
	x = x_start;
	t = 0.8 * (n - 1);
	output.open("wykres_t1_6.txt", ios::out);
	for(int i = 0; i < m; i+=5){
		output << setprecision(8) << x << " " << analytical[t][i] << " ";
		output << thomas_l[t][i] << " " << thomas_cn[t][i] << " ";
		output << gauss_l[t][i] << " " << gauss_cn[t][i] << endl;
		x += h;
	}
	output.close();

	// t = 2s
	x = x_start;
	t = n - 1;
	output.open("wykres_t2.txt", ios::out);
	for(int i = 0; i < m; i+=5){
		output << setprecision(8) << x << " " << analytical[t][i] << " ";
		output << thomas_l[t][i] << " " << thomas_cn[t][i] << " ";
		output << gauss_l[t][i] << " " << gauss_cn[t][i] << endl;
		x += h;
	}
	output.close();
}


/**
 * Funkcja zapisująca do pliku wartość masymalnych błędów bezwględnych w funkcji czasu
 * @param analytical macierzy wyników analitycznych
 * @param thomas_l macierzy wyników numerycznych dla alg Thomasa i dyskretyzacji Laasonen
 * @param thomas_cn macierz wyników numerycznych dla alg Thomasa i dyskretyzacji Cranka-Nicolson
 * @param gauss_l macierzy wyników numerycznych dla metody Gaussa-Seidela i dyskretyzacji Laasonen
 * @param gauss_cn macierzy wyników numerycznych dla metody Gaussa-Seidela i dyskretyzacji Cranka-Nicolson
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn macierzy
 * @param delta_t krok czasowy delta_t
 */
void time_error(double **analytical, double **thomas_l, double **thomas_cn, double **gauss_l, double **gauss_cn, int n, int m, double delta_t){
	//wyliczamy kolejne wartości t
	double *step_t = get_steps_t(delta_t, n, t_0);

	//Wektor błędów dla Thomasa + Laasonen
	double **error_matrix = get_error_matrix(analytical, thomas_l, n, m);
	double *error_vector_tl = get_error_vector(error_matrix, n, m);
	delete_matrix(error_matrix, n);

	///Wektor błędów dla Thomasa + CN
	error_matrix = get_error_matrix(analytical, thomas_cn, n, m);
	double *error_vector_tcn = get_error_vector(error_matrix, n, m);
	delete_matrix(error_matrix, n);

	//Wektor błędów dla Gaussa + Laasonen
	error_matrix = get_error_matrix(analytical, gauss_l, n, m);
	double *error_vector_gl = get_error_vector(error_matrix, n, m);
	delete_matrix(error_matrix, n);

	//Wektor błędów dla Gaussa + CN
	error_matrix = get_error_matrix(analytical, gauss_cn, n, m);
	double *error_vector_gcn = get_error_vector(error_matrix, n, m);
	delete_matrix(error_matrix, n);

	fstream output;
	//W pliku zapisujemy wartości t oraz błędów
	output.open("time_errors.txt", ios::out);
	for(int i = 0; i < n; i++){
		output << setprecision(8) << step_t[i] << " " << error_vector_tl[i] << " ";
		output << error_vector_tcn[i] << " " << error_vector_gl[i] << " ";
		output << error_vector_gcn[i] << endl;
	}
	output.close();
	
	//usuwamy zaalokowane wektory
	delete [] step_t;
	delete [] error_vector_tl;
	delete [] error_vector_tcn;
	delete [] error_vector_gl;
	delete [] error_vector_gcn;
}


/**
 * Funkcja wykonująca potrzebne obliczenia
 */
void calculate_all(){
	fstream output;

	output.open("max_error_t.txt", ios::out);
	//Dla każdego h liczymy wartości numeryczne
	for(int m = 100; m <= 1000; m += 100){
		double h = (x_end - x_start) / (m - 1);
		double delta_t = get_delta_t(LAMBDA, h, D);
		int n = t_max / delta_t + 1;
		cout << n << " x " << m  << " " << log10(h) << endl;
 
		double **analytical = get_analytical_solution(n, m, h, delta_t);
		double **thomas_l = get_Laasonen_Thomas(n, m);
		double **thomas_cn = get_Crank_Nicolson_Thomas(n, m);
		double **gauss_l = get_Laasonen_Gauss_Seidel(n, m);
		double **gauss_cn = get_Crank_Nicolson_Gauss_Seidel(n, m);

		//Wyliczamy wektor błędów bezwzględnych dla danej metody
		double *thomas_l_err = error_vector(analytical[n - 1], thomas_l[n - 1], m);
		double *thomas_cn_err = error_vector(analytical[n - 1], thomas_cn[n - 1], m);
		double *gauss_l_err = error_vector(analytical[n - 1], gauss_l[n - 1], m);
		double *gauss_cn_err = error_vector(analytical[n - 1], gauss_cn[n - 1], m);

		//W pliku zapisujemy logarytm dziesiętny z błędu maskymalnego danego wektora oraz h
		output << setprecision(15) << log10(h) << " " << get_max_error(thomas_l_err, m) << " ";
		output << get_max_error(thomas_cn_err, m) << " ";
		output << get_max_error(gauss_l_err, m) << " " << get_max_error(gauss_cn_err, m) << endl;

		//Dla optymalnego h, to znaczy przy 800 węzłach zapisujemy wyniki dla punktu 2 i 3
		//zapisujemy wyniki analityczne, numeryczne i błędy
		if(m == 800){
			draw_plot_time(analytical, thomas_l, thomas_cn, gauss_l, gauss_cn, n, m, 5 * h);
			time_error(analytical, thomas_l, thomas_cn, gauss_l, gauss_cn, n, m, delta_t);
			write_results(analytical[n - 1], thomas_l[n - 1], thomas_l_err, m, h, x_start, "thomas_laasonen.txt");
			write_results(analytical[n - 1], thomas_cn[n - 1], thomas_cn_err, m, h, x_start, "thomas_crank_nicolson.txt");
			write_results(analytical[n - 1], gauss_l[n - 1], gauss_l_err, m, h, x_start, "gauss_seidel_laasonen.txt");
			write_results(analytical[n - 1], gauss_cn[n - 1], gauss_cn_err, m, h, x_start, "gauss_seidel_crank_nicolson.txt");
		}

		//Usuwamy zaalokowane zasoby
		delete [] thomas_l_err;
		delete [] thomas_cn_err;
		delete [] gauss_l_err;
		delete [] gauss_cn_err;

		clean_up(analytical, thomas_l, thomas_cn, gauss_l, gauss_cn, n);
	}

	output.close();
}


/**
 * Funkcja main
 * @return 0 jeśli zadziała 
 */
int main(){
	calculate_all();

	return 0;
}
