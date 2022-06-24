#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#define TOLF 1e-15 //zadana tolerancja residuum w metodzie Gaussa-Seidela
#define TOLX 1e-15 //zadana tolerancja błędu

using namespace std;

/**
 * Funkcja zwracająca macierz z warunkami początkowymi i brzegowymi
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn macierzy
 * @return double** zaalokowana macierz
 */
double **get_matrix_conditions(int n, int m);


/**
 * Funkcja zwaracająca obliczoną wartość kroku czasowego
 * @param lambda parametr lambda metod pośrednich
 * @param h krok przestrzenny
 * @param D współczynnik transportu ciepła
 * @return double obliczony krok delta_t
 */
double get_delta_t(const double lambda, const double h, const double D);


/**
 * Funkcja pomocniczna do alokacji macierzy
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn
 * @return double** zaalokowana macierz n x m
 */
double **get_matrix(int n, int m);


/**
 * Funkcja pomocnicza alokująca wektor
 * @param m rozmiar wektora
 * @return double* zaalokowany wektor
 */
double *get_vector(int m);


/**
 * Funkcja obliczająca macierz błędów bezwzględnych
 * @param analytical macierz wyników analitycznych
 * @param result macierz wyników obliczonych numerycznie
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn macierzy
 * @return double** macierz błędów względnych
 */
double **get_error_matrix(double **analytical, double **result, const int n, const int m);


/**
 * Funkcja obliczająca wektor błędów maksymalnych w danym wierszu macierzy
 * @param err macierz błędów
 * @param n liczba wierszy macierzy
 * @param m liczba kolumn macierzy
 * @return double* wektor błędu maksymalnego
 */
double *get_error_vector(double **err, const int n, const int m);


/**
 * Funkcja licząca wektor błędu dla poziomu czasowego t = t_max
 * @param a wektor rozwiązań analitycznych
 * @param x wektor rozwiązań numerycznych
 * @param m rozmiar wektora
 * @return double* wektor błędów wzglęnych
 */
double *error_vector(double *a, double *x, int m);


/**
 * Funkcja zwracająca logarytm dziesiętny z maksymalnego błędu w wektorze
 * @param v wektor błędu
 * @param m rozmiar wektora
 * @return double logarytm dziesiętny ze znalezionego maksimum
 */
double get_max_error(double *v, const int m);


/**
 * Funkcja obliczająca kolejne wartości czasu t
 * @param delta_t krok czasowy
 * @param n rozmiar wektora
 * @param t_0 punkt początkowy
 * @return double* wektor wartości czasu
 */
double *get_steps_t(double delta_t, int n, double t_0);


/**
 * Funkcja usuwająca zaalokowane macierze
 * @param results macierz wyników
 * @param a macierz do usunięcia
 * @param b macierz do usunięcia 
 * @param c macierz do usunięcia
 * @param d macierz do usunięcia
 * @param n liczba wierszy
 */
void clean_up(double **results, double **a, double **b, double **c, double **d, int n);


/**
 * Funkcja usuwająca macierz
 * @param matrix usuwana macierz
 * @param n liczba wierszy
 */
void delete_matrix(double **matrix, int n);


/**
 * Funkcja zapisująca do pliku wartości analityczne, numeryczne i odpowiadające
 * im błędy bezwzględne dla optymalnego h w chwili t_max
 * @param analitical wyniki analityczne
 * @param numerical wyniki numeryczne
 * @param err błędy bezwzględne
 * @param m rozmiar wektorów
 * @param h krok przestrzenny
 * @param filename nazwa pliku
 */
void write_results(double *analitical, double *numerical, double *err, int m, double h, double x_start, const char* filename);


/**
 * Funkcja sprawdzająca kryterium dokładności wyznaczenia x_n
 * @param x_n obliczony wektor x_n
 * @param x_n1 wektor z poprzeniego przybliżenia
 * @param m rozmiar wektorów
 * @return true jeśli wektor mieści się w zadanej tolerancji błędu
 * @return false jeśli wektor nie spełnia kryterium
 */
bool check_en(double *x_n, double* x_n1, int m);


/**
 * Funkcja sprawdzająca kryterium wiarygodności
 * @param l wektor wartości pod przekątną główną
 * @param d wektor wartości na przekątnej głównej
 * @param u wektor wartości nad przekątną główną
 * @param x wektor wartości przybliżonych
 * @param b wektor wyniku
 * @param m rozmiar wektorów
 * @return true jeśli wynik mieści się w zadanej tolerancji
 * @return false jesli wynik nie mieści się w zadanej tolerancji
 */
bool check_residue(double* l, double *d, double *u, double *x, double* b, int m);