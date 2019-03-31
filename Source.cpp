#include "Element.h"

void load(float&, float&, float&, float&, float&, float&, float&, int&, int&, float&, float&, float&);
void grid(float, float, float, float, float, float, float, int, int, float, float, float, Node*, Element*);
void tab0(float**&, int);
void vec0(float*&, int);
void print(float**, int);
void agregationH(float**&, Element*);
void agregationC(float**&, Element*);
void agregationHBC(float**&, Element*);
void agregationP(float*&, Element*);
bool ludist(int, float**);
bool lusolve(int, int, float**, float**);
void inverse(float**, float**&, int);
void mnoz(float *&, float**, float*, int);
void vecP(float*&, int, float, float**, float**, Element*);
void MIN_MAX(float*, int, float&, float&);
const double eps = 1e-12;

int main()
{
	int N_H, N_B;
	float initial_temperature, simulation_time, simulation_step_time, To, alfa, H, B, c, k, ro;
	load(initial_temperature, simulation_time, simulation_step_time, To, alfa, H, B, N_H, N_B, c, k, ro);
	int bok = N_H*N_B;
	int quantity = (N_H - 1)*(N_B - 1);
	float **HG;
	float **CG;
	float **HBCG;
	float **HHBC;
	float **Z;
	float **Z2;
	float **CDT;
	float *P;
	vec0(P, bok);
	tab0(HG, bok);
	tab0(CG, bok);
	tab0(HBCG, bok);
	tab0(HHBC, bok);
	tab0(Z, bok);
	tab0(CDT, bok);
	tab0(Z2, bok);
	Element *element = new Element[quantity];
	Node *node = new Node[bok];
	grid(initial_temperature, simulation_time, simulation_step_time, To, alfa, H, B, N_H, N_B, c, k, ro, node, element);
	agregationH(HG, element);
	agregationC(CG, element);
	agregationHBC(HBCG, element);
	agregationP(P, element);
	for (int i = 0; i < bok; i++)
		for (int j = 0; j < bok; j++)
			HHBC[i][j] = HBCG[i][j] + HG[i][j] + (CG[i][j] / simulation_step_time);
	for (int i = 0; i < bok; i++)
		for (int j = 0; j < bok; j++)
			Z[i][j] = HG[i][j] + HBCG[i][j];
	for (int i = 0; i < bok; i++)
		for (int j = 0; j < bok; j++)
			CDT[i][j] = CG[i][j] / simulation_step_time;
	cout << "__________________________________________________MATRIX [C]__________________________________________________" << endl;
	print(CG, bok);
	cout << "__________________________________________________MATRIX [H]__________________________________________________" << endl;
	print(HG, bok);
	cout << "______________________________________________MATRIX ([H]+[C]/dT)_____________________________________________" << endl;
	print(HHBC, bok);
	vecP(P, 10, initial_temperature, HHBC, CDT, element);
	system("pause");
}
void load(float &initial_temperature, float &simulation_time, float &simulation_step_time, float &To, float &alfa, float &H, float &B, int &N_H, int &N_B, float &c, float &k, float &ro)
{
	string tekst[12];
	int i = 0;
	fstream plik;
	plik.open("plik.txt", ios::in);
	if (plik.good() == true)
	{
		while (!plik.eof())
		{
			getline(plik, tekst[i]);
			i++;
		}
	}
	else
		cout << "ERROR" << endl;
	plik.close();
	initial_temperature = atof(tekst[0].c_str());
	simulation_time = atof(tekst[1].c_str());
	simulation_step_time = atof(tekst[2].c_str());
	To = atof(tekst[3].c_str());
	alfa = atof(tekst[4].c_str());
	H = atof(tekst[5].c_str());
	B = atof(tekst[6].c_str());
	N_H = atoi(tekst[7].c_str());
	N_B = atoi(tekst[8].c_str());
	c = atof(tekst[9].c_str());
	k = atof(tekst[10].c_str());
	ro = atof(tekst[11].c_str());
}
void grid(float initial_temperature, float simulation_time, float simulation_step_time, float To, float alfa, float H, float B, int N_H, int N_B, float c, float k, float ro, Node *node, Element *element)
{
	if (N_B > 0 && N_H >0) {


		double delta_y = H / (N_H - 1);
		double delta_x = B / (N_B - 1);
		int kolumna = 1;
		int przesuniecie = 0;
		for (int i = 0; i<N_H*N_B; i++) 
		{
			if ((i*1.0) / N_H >= 1 && i%N_H == 0)
				kolumna++;

			node[i].set(delta_x*(kolumna - 1), delta_y*(i%N_H));
		}

		for (int i = 0; i<(N_H - 1)*(N_B - 1); i++) {
			if ((i + przesuniecie) % N_H == (N_H - 1))
				++przesuniecie;
			element[i].set(i + przesuniecie, node, initial_temperature, simulation_time, simulation_step_time, To, alfa, H, B, N_H, N_B, c, k, ro);
		}
	}
	else {
		cout << "DANE NIEPOPRAWNIE WCZYTANE Z PLIKU!!" << endl;
	}
}
void tab0(float **&tab, int x)
{
	tab = new float*[x];
	for (int i = 0; i < x; i++)
	{
		tab[i] = new float[x];
	}
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x; j++)
		{
			tab[i][j] = 0;
		}
	}
}
void vec0(float *&tab, int x)
{
	tab = new float[x];
	for (int i = 0; i < x; i++)
		tab[i] = 0;
}
void print(float **tab, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x; j++)
		{
			cout << tab[i][j] << "  ";
		}
		cout << endl;
	}
}
void agregationH(float **&tab, Element *element)
{
	for (int a = 0; a < (element[a].get_nh() - 1)*(element[a].get_nb() - 1); a++)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				tab[element[a].get_id(i)][element[a].get_id(j)] += element[a].get_matrixh(i, j);
}
void agregationC(float **&tab, Element *element)
{
	for (int a = 0; a < (element[a].get_nh() - 1)*(element[a].get_nb() - 1); a++)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				tab[element[a].get_id(i)][element[a].get_id(j)] += element[a].get_matrixc(i, j);
}
void agregationHBC(float **&tab, Element *element)
{
	for (int a = 0; a < (element[a].get_nh() - 1)*(element[a].get_nb() - 1); a++)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				tab[element[a].get_id(i)][element[a].get_id(j)] += element[a].get_matrixhbc(i, j);
}
void agregationP(float *&tab, Element *element)
{
	for (int a = 0; a < (element[a].get_nh()-1)*(element[a].get_nb() - 1); a++)
		for (int i = 0; i < 4; i++)
			tab[element[a].get_id(i)] = element[a].get_p(i);
}
bool ludist(int n, float **A)
{
	int i, j, k;

	for (k = 0; k < n - 1; k++)
	{
		if (fabs(A[k][k]) < eps) return false;

		for (i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
	}

	return true;
}
bool lusolve(int k, int n, float **A, float **X)
{
	int    i, j;
	double s;

	for (i = 1; i < n; i++)
	{
		s = 0;

		for (j = 0; j < i; j++) s += A[i][j] * X[j][k];

		X[i][k] -= s;
	}

	if (fabs(A[n - 1][n - 1]) < eps) return false;

	X[n - 1][k] /= A[n - 1][n - 1];

	for (i = n - 2; i >= 0; i--)
	{
		s = 0;

		for (j = i + 1; j < n; j++) s += A[i][j] * X[j][k];

		if (fabs(A[i][i]) < eps) return false;

		X[i][k] = (X[i][k] - s) / A[i][i];
	}

	return true;
}
void inverse(float **A, float **&X, int n)
{
	X = new float*[n];
	for (int a = 0; a < n; a++)
		X[a] = new float[n];
	int i, j;
	bool ok;
	if (ludist(n, A))
	{  
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++) X[i][j] = 0;
			X[i][i] = 1;
		}

		ok = true;
		for (i = 0; i < n; i++)
			if (!lusolve(i, n, A, X))
			{
				ok = false;
				break;
			}
	}
	else ok = false;
}
void mnoz(float *&T, float **A, float *B, int n)
{
	for (int i = 0; i < 16; i++)
		T[i] = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			T[i] += (A[i][j] * B[j]);
}
void vecP(float *&P, int iter, float initial_temperature, float **HHBC, float **CDT, Element *element)
{
	float MAX = 0, MIN = 0;
	float **X;
	inverse(HHBC, X, 16);
	float *T = new float[16];
	for (int i = 0; i < 16; i++)
		T[i] = initial_temperature;
	float *T2 = new float[16];
	for (int i = 0; i < 16; i++)
		T2[i] = 0;

	for (int i = 0; i < iter; i++)
	{
		if (i == 0)
		{
			

			mnoz(T2, CDT, T, 16);
			for (int i = 0; i < 16; i++)
				P[i] = P[i] + T2[i];
			cout << "__________________________________________Vector ([{P}+{[C]/dT}*{T0})_________________________________________" << endl;
			cout << "__________________________________________________iteration "<<i<<"_________________________________________________" << endl;
			for (int i = 0; i < 16; i++)
				cout << P[i] << "  ";
			cout << endl << endl;
			mnoz(T, X, P, 16);
			MIN_MAX(T, 16, MIN, MAX);
			cout << "STEP: " << (i + 1) * 50 << ",   TAMP MAX: "<< MIN << ",   TEMP MIN:" << MAX << endl;
		}
		else
		{
			
			mnoz(T2, CDT, T, 16);
			for (int i = 0; i < 16; i++)
				P[i] = 0;
			for (int a = 0; a < 9; a++)
				for (int i = 0; i < 4; i++)
					P[element[a].get_id(i)] = element[a].get_p(i);
			for (int i = 0; i < 16; i++)
				P[i] = P[i] + T2[i];
			cout << "__________________________________________Vector ([{P}+{[C]/dT}*{T0})_________________________________________" << endl;
			cout << "__________________________________________________iteration " << i << "_________________________________________________" << endl;
			for (int i = 0; i < 16; i++)
				cout << P[i] << "  ";
			cout << endl << endl;
			mnoz(T, X, P, 16);
			MIN_MAX(T, 16, MIN, MAX);
			cout <<"STEP: "<<(i+1)*50<< ",   TEMP MIN: " << MIN << ",   TEMP MAX: " << MAX << endl;
		}
	}
}
void MIN_MAX(float *tab, int x, float &MIN, float &MAX)
{
	MAX = tab[0];
	MIN = tab[0];
	for (int i = 0; i < x; i++)
	{
		if (tab[i] > MAX)
			MAX = tab[i];
		if (tab[i] < MIN)
			MIN = tab[i];
	}
}