#include "Element.h"
Element::Element()
{

}
void Element::set(int x, Node *node, float Initial_temperature, float Simulation_time, float Simulation_step_time, float to, float Alfa, float h, float b, int n_H, int n_B, float C, float K, float Ro)
{
	ID = new int[4];
	ID[0] = x;
	ID[1] = x + n_H;
	ID[2] = x + n_H + 1;
	ID[3] = x + 1;
	node1 = node[x];
	node2 = node[x + n_H];
	node3 = node[x + n_H + 1];
	node4 = node[x + 1];
	initial_temperature = Initial_temperature;
	simulation_time = Simulation_step_time;
	simulation_time = Simulation_time;
	To = to;
	alfa = Alfa;
	H = h;
	B = b;
	N_H = n_H;
	N_B = n_B;
	c = C;
	k = K;
	ro = Ro;
	if (ID[0]<(N_H - 1))
		flagal = 1;
	if (ID[1] > ((N_H*N_B) - N_H - 1))
		flagap = 1;
	if (ID[0] % N_H == 0)
		flagad = 1;
	if (ID[3] % N_H == (N_H - 1))
		flagag = 1;
	fun();
	matrix();
	matrixc();
	matrixhbc();
	p();

}
void Element::fun()
{
	float *ksi = new float[4];
	ksi[0] = (-1) / sqrt(3);
	ksi[1] = (-1)*ksi[0];
	ksi[2] = ksi[1];
	ksi[3] = (-1)*ksi[2];
	float *eta = new float[4];
	eta[0] = ksi[0];
	eta[1] = eta[0];
	eta[2] = (-1)*eta[1];
	eta[3] = eta[2];
	float *N1 = new float[4];
	float *N2 = new float[4];
	float *N3 = new float[4];
	float *N4 = new float[4];
	for (int i = 0; i < 4; i++)
	{
		N1[i] = (0.25)*(1 - ksi[i])*(1 - eta[i]);
		N2[i] = (0.25)*(1 + ksi[i])*(1 - eta[i]);
		N3[i] = (0.25)*(1 + ksi[i])*(1 + eta[i]);
		N4[i] = (0.25)*(1 - ksi[i])*(1 + eta[i]);
	}
	float *XP = new float[4];
	float *YP = new float[4];
	for (int i = 0; i < 4; i++)
	{
		XP[i] = (N1[i] * node1.get_x()) + (N2[i] * node2.get_x()) + (N3[i] * node3.get_x()) + (N4[i] * node4.get_x());
		YP[i] = (N1[i] * node1.get_y()) + (N2[i] * node2.get_y()) + (N3[i] * node3.get_y()) + (N4[i] * node4.get_y());
	}
	float tab1[4][4];
	float tab2[4][4];
	for (int i = 0; i < 4; i++)
	{
		tab1[0][i] = (-0.25)*(1 - eta[i]);
		tab1[1][i] = (0.25)*(1 - eta[i]);
		tab1[2][i] = (0.25)*(1 + eta[i]);
		tab1[3][i] = (-0.25)*(1 + eta[i]);
		tab2[0][i] = (-0.25)*(1 - ksi[i]);
		tab2[1][i] = (-0.25)*(1 + ksi[i]);
		tab2[2][i] = (0.25)*(1 + ksi[i]);
		tab2[3][i] = (0.25)*(1 - ksi[i]);

	}
	jakobian1 = new float*[4];
	for (int i = 0; i < 4; i++)
		jakobian1[i] = new float[4];
	for (int i = 0; i < 4; i++)
	{
		jakobian1[0][i] = (tab1[0][i] * node1.get_x()) + (tab1[1][i] * node2.get_x()) + (tab1[2][i] * node3.get_x()) + (tab1[3][i] * node4.get_x());
		jakobian1[1][i] = (tab1[0][i] * node1.get_y()) + (tab1[1][i] * node2.get_y()) + (tab1[2][i] * node3.get_y()) + (tab1[3][i] * node4.get_y());
		jakobian1[2][i] = (tab2[0][i] * node1.get_x()) + (tab2[1][i] * node2.get_x()) + (tab2[2][i] * node3.get_x()) + (tab2[3][i] * node4.get_x());
		jakobian1[3][i] = (tab2[0][i] * node1.get_y()) + (tab2[1][i] * node2.get_y()) + (tab2[2][i] * node3.get_y()) + (tab2[3][i] * node4.get_y());
	}
	detJ = new float[4];
	for (int i = 0; i < 4; i++)
	{
		detJ[i] = (jakobian1[0][i] * jakobian1[3][i]) - (jakobian1[1][i] * jakobian1[2][i]);
	}
	jakobian = new float*[4];
	for (int i = 0; i < 4; i++)
		jakobian[i] = new float[4];
	for (int i = 0; i < 4; i++)
	{
		jakobian[0][i] = jakobian1[3][i] / detJ[i];
		jakobian[1][i] = (-1)*jakobian1[1][i] / detJ[i];
		jakobian[2][i] = jakobian1[2][i] / detJ[i];
		jakobian[3][i] = jakobian1[0][i] / detJ[i];

	}
}
void Element::matrix()
{
	float *ksi = new float[4];
	ksi[0] = (-1) / sqrt(3);
	ksi[1] = (-1)*ksi[0];
	ksi[2] = ksi[1];
	ksi[3] = (-1)*ksi[2];
	float *eta = new float[4];
	eta[0] = ksi[0];
	eta[1] = eta[0];
	eta[2] = (-1)*eta[1];
	eta[3] = eta[2];
	float pierw = 1 / sqrt(3);
	float tab1[4][4];
	float tab2[4][4];
	float **tab3 = new float *[4];
	for (int i = 0; i<4; i++)
		tab3[i] = new float[4];
	float **tab4 = new float *[4];
	for (int i = 0; i<4; i++)
		tab4[i] = new float[4];
	for (int i = 0; i < 4; i++)
	{
		tab1[0][i] = (-0.25)*(1 - eta[i]);
		tab1[1][i] = (0.25)*(1 - eta[i]);
		tab1[2][i] = (0.25)*(1 + eta[i]);
		tab1[3][i] = (-0.25)*(1 + eta[i]);
		tab2[0][i] = (-0.25)*(1 - ksi[i]);
		tab2[1][i] = (-0.25)*(1 + ksi[i]);
		tab2[2][i] = (0.25)*(1 + ksi[i]);
		tab2[3][i] = (0.25)*(1 - ksi[i]);
	}
	for (int i = 0; i < 4; i++)
	{
		tab3[0][i] = (jakobian[0][0] * tab1[i][0]) + (jakobian[1][0] * tab2[i][0]);
		tab3[1][i] = (jakobian[0][1] * tab1[i][1]) + (jakobian[1][1] * tab2[i][1]);
		tab3[2][i] = (jakobian[0][2] * tab1[i][2]) + (jakobian[1][2] * tab2[i][2]);
		tab3[3][i] = (jakobian[0][3] * tab1[i][3]) + (jakobian[1][3] * tab2[i][3]);
		tab4[0][i] = (jakobian[2][0] * tab1[i][0]) + (jakobian[3][0] * tab2[i][0]);
		tab4[1][i] = (jakobian[2][1] * tab1[i][1]) + (jakobian[3][1] * tab2[i][1]);
		tab4[2][i] = (jakobian[2][2] * tab1[i][2]) + (jakobian[3][2] * tab2[i][2]);
		tab4[3][i] = (jakobian[2][3] * tab1[i][3]) + (jakobian[3][3] * tab2[i][3]);
	}
	float ***tab5 = new float **[8];
	for (int i = 0; i < 8; i++)
	{
		tab5[i] = new float *[4];
		for (int j = 0; j < 4; j++)
		{
			tab5[i][j] = new float[4];
		}
	}
	fun1(tab5, tab3, 0, 0);
	fun1(tab5, tab3, 1, 1);
	fun1(tab5, tab3, 2, 2);
	fun1(tab5, tab3, 3, 3);
	fun1(tab5, tab4, 4, 0);
	fun1(tab5, tab4, 5, 1);
	fun1(tab5, tab4, 6, 2);
	fun1(tab5, tab4, 7, 3);
	det(tab5);
	float ***tab6 = new float **[4];
	for (int i = 0; i < 8; i++)
	{
		tab6[i] = new float *[4];
		for (int j = 0; j < 4; j++)
		{
			tab6[i][j] = new float[4];
		}
	}
	tab6 = fun2(tab5);
	matrixH = new float*[4];
	for (int i = 0; i < 4; i++)
		matrixH[i] = new float[4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrixH[i][j] = tab6[0][i][j] + tab6[1][i][j] + tab6[2][i][j] + tab6[3][i][j];
		}
	}
}
void Element::matrixc()
{
	float *ksi = new float[4];
	ksi[0] = (-1) / sqrt(3);
	ksi[1] = (-1)*ksi[0];
	ksi[2] = ksi[1];
	ksi[3] = (-1)*ksi[2];
	float *eta = new float[4];
	eta[0] = ksi[0];
	eta[1] = eta[0];
	eta[2] = (-1)*eta[1];
	eta[3] = eta[2];
	float *N1 = new float[4];
	float *N2 = new float[4];
	float *N3 = new float[4];
	float *N4 = new float[4];
	for (int i = 0; i < 4; i++)
	{
		N1[i] = (0.25)*(1 - ksi[i])*(1 - eta[i]);
		N2[i] = (0.25)*(1 + ksi[i])*(1 - eta[i]);
		N3[i] = (0.25)*(1 + ksi[i])*(1 + eta[i]);
		N4[i] = (0.25)*(1 - ksi[i])*(1 + eta[i]);
	}
	float ***tab = new float **[4];
	for (int i = 0; i < 8; i++)
	{
		tab[i] = new float *[4];
		for (int j = 0; j < 4; j++)
		{
			tab[i][j] = new float[4];
		}
	}
	float cro = mnoz();
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab[0][i][j] = N1[j] * N1[i] * detJ[0] * cro;
			tab[1][i][j] = N2[j] * N2[i] * detJ[1] * cro;
			tab[2][i][j] = N3[j] * N3[i] * detJ[2] * cro;
			tab[3][i][j] = N4[j] * N4[i] * detJ[3] * cro;
		}
	}
	matrixC = new float*[4];
	for (int i = 0; i < 4; i++)
		matrixC[i] = new float[4];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrixC[i][j] = tab[0][i][j] + tab[1][i][j] + tab[2][i][j] + tab[3][i][j];
		}
	}
}
void Element::matrixhbc()
{
	float *ksi = new float[2];
	ksi[0] = (-1) / sqrt(3);
	ksi[1] = (-1)*ksi[0];
	float *N1 = new float[2];
	float *N2 = new float[2];
	for (int i = 0; i < 2; i++)
	{
		N1[i] = 0.5*(1 - ksi[i]);
		N2[i] = 0.5*(1 + ksi[i]);
	}
	float bok = H / (N_H-1);
	float Detj = bok / 2;
	float tab1[2][2];
	float tab2[2][2];
	float tab3[2][2];
	float pow1[2][2];
	float pow2[2][2];
	float pow3[2][2];
	float pow4[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			tab1[i][j] = N1[i] * N1[j] * alfa;
			tab2[i][j] = N2[i] * N2[j] * alfa;
			tab3[i][j] = (tab1[i][j] + tab2[i][j])*Detj;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			pow1[i][j] = tab3[i][j] * int(flagad);
			pow2[i][j] = tab3[i][j] * int(flagap);
			pow3[i][j] = tab3[i][j] * int(flagag);
			pow4[i][j] = tab3[i][j] * int(flagal);
		}
	}
	matrixHBC = new float*[4];
	for (int i = 0; i < 4; i++)
		matrixHBC[i] = new float[4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			matrixHBC[i][j] = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i < 2 && j < 2)
				matrixHBC[i][j] += pow1[i][j];
			if (i > 0 && i < 3 && j < 3 && j>0)
				matrixHBC[i][j] += pow2[i - 1][j - 1];
			if (i > 1 && j>1)
				matrixHBC[i][j] += pow3[i - 2][j - 2];
			if (i == 0 && j == 0)
				matrixHBC[i][j] += pow4[i][j];
			if (i == 3 && j == 3)
				matrixHBC[i][j] += pow4[i - 2][j - 2];
			if (i == 3 && j == 0)
				matrixHBC[i][j] += pow4[i - 2][j];
			if (i == 0 && j == 3)
				matrixHBC[i][j] += pow4[i][j - 2];
		}
	}
	/*matrixHBC[0][0] = pow1[0][0] + pow4[0][0];
	matrixHBC[0][1] = pow1[0][1];
	matrixHBC[0][2] = 0;
	matrixHBC[0][3] = pow4[0][1];
	matrixHBC[1][0] = pow1[1][0];
	matrixHBC[1][1] = pow1[1][1] + pow2[0][0];
	matrixHBC[1][2] = pow2[0][1];
	matrixHBC[1][3] = 0;
	matrixHBC[2][0] = 0;
	matrixHBC[2][1] = pow2[1][0];
	matrixHBC[2][2] = pow2[1][1] + pow3[0][0];
	matrixHBC[2][3] = pow3[0][1];
	matrixHBC[3][0] = pow4[1][0];
	matrixHBC[3][1] = 0;
	matrixHBC[3][2] = pow3[1][0];
	matrixHBC[3][3] = pow3[1][1] + pow4[1][1];*/


}
void Element::p()
{
	float *ksi = new float[2];
	ksi[0] = (-1) / sqrt(3);
	ksi[1] = (-1)*ksi[0];
	float *N1 = new float[2];
	float *N2 = new float[2];
	for (int i = 0; i < 2; i++)
	{
		N1[i] = 0.5*(1 - ksi[i]);
		N2[i] = 0.5*(1 + ksi[i]);
	}
	float deltax = H / (N_H-1);
	float Detj = deltax / 2;
	float sum;
	sum = (N1[0] * alfa*Detj*To) + (N2[0] * alfa*Detj*To) + (N1[1] * alfa*Detj*To) + (N2[1] * alfa*Detj*To);
	P = new float[4];
	bool F[4];
	for (int i = 0; i < 4; i++)
		F[i] = 0;
	if (flagal == 1 || flagad == 1)
		F[0] = 1;
	if (flagad == 1 || flagap == 1)
		F[1] = 1;
	if (flagap == 1 || flagag == 1)
		F[2] = 1;
	if (flagag == 1 || flagal == 1)
		F[3] = 1;

	P[0] = int(F[0])*sum;
	P[1] = int(F[1])*sum;
	P[2] = int(F[2])*sum;
	P[3] = int(F[3])*sum;

}
void Element::fun1(float ***&tab, float **&tab2, int x, int y)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab[x][i][j] = tab2[y][j] * tab2[y][i];
		}
	}
}
void Element::det(float ***&tab)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab[0][i][j] = tab[0][i][j] * detJ[0];
			tab[1][i][j] = tab[1][i][j] * detJ[1];
			tab[2][i][j] = tab[2][i][j] * detJ[2];
			tab[3][i][j] = tab[3][i][j] * detJ[3];
			tab[4][i][j] = tab[4][i][j] * detJ[0];
			tab[5][i][j] = tab[5][i][j] * detJ[1];
			tab[6][i][j] = tab[6][i][j] * detJ[2];
			tab[7][i][j] = tab[7][i][j] * detJ[3];
		}
	}
}
float ***Element::fun2(float ***tab)
{
	float ***tab1 = new float **[4];
	for (int i = 0; i < 8; i++)
	{
		tab1[i] = new float *[4];
		for (int j = 0; j < 4; j++)
		{
			tab1[i][j] = new float[4];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab1[0][i][j] = (tab[0][i][j] + tab[4][i][j])*k;
			tab1[1][i][j] = (tab[1][i][j] + tab[5][i][j])*k;
			tab1[2][i][j] = (tab[2][i][j] + tab[6][i][j])*k;
			tab1[3][i][j] = (tab[3][i][j] + tab[7][i][j])*k;
		}
	}
	return tab1;
}
float Element::mnoz()
{
	int x = c;
	int y = ro;
	return x*y;
}
void Element::set_detJ(float *detj)
{
detJ = detj;
}
void Element::set_Jakobian(float **Jakobian)
{
jakobian = Jakobian;
}
void Element::set_Jakobian1(float **Jakobian1)
{
jakobian1 = Jakobian1;
}
void Element::set_matrixH(float **matrixh)
{
	matrixH = matrixh;
}
void Element::set_matrixC(float **matrixc)
{
	matrixC = matrixc;
}
void Element::set_matrixHBC(float **matrixhbc)
{
	matrixHBC = matrixhbc;
}
void Element::set_flagal(bool flaga)
{
	flagal = flaga;
}
void Element::set_flagap(bool flaga)
{
	flagap = flaga;
}
void Element::set_flagag(bool flaga)
{
	flagag = flaga;
}
void Element::set_flagad(bool flaga)
{
	flagad = flaga;
}
int Element::get_id(int i)
{
	return ID[i];
}
float Element::get_matrixh(int i, int j)
{
	return matrixH[i][j];
}
float Element::get_matrixc(int i, int j)
{
	return matrixC[i][j];
}
float Element::get_matrixhbc(int i, int j)
{
	return matrixHBC[i][j];
}
int Element::get_p(int i)
{
	return P[i];
}
float Element::get_nh()
{
	return N_H;
}
float Element::get_nb()
{
	return N_B;
}
Element::~Element()
{
	for (int i = 0; i < N_H; i++)
	{
		delete[] jakobian[i];
		delete[] jakobian1[i];
		delete[] matrixH[i];
		delete[] matrixC[i];
		delete[] matrixHBC[i];
	}
	delete[] ID;
	delete[] jakobian;
	delete[] jakobian1;
	delete[] matrixH;
	delete[] matrixC;
	delete[] matrixHBC;
	delete[] detJ;
	delete[] P;
		
}
