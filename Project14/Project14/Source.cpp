/***
Лабораторная работа 2 по дисциплине "Модели решения задач в интеллектуальных системах"
Выполнено студентами БГУИР группы 821701
Ясинский Никита Максимович
Трипутько Роман Викторович
Модель решателя задач с ОКМД архитектурой
04.05.2020
Класс ОКМД, реализующий модель архитектуры ОКМД
*/


#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;


class OKMD {

private:
	vector<vector<double>> matrA;
	vector<vector<double>> matrB;
	vector<vector<double>> matrE;
	vector<vector<double>> matrG;
	vector<vector<double>> matrC;

	vector<vector<vector<double>>> matrD;
	vector<vector<vector<double>>> matrF;

	const double MIN = -1;
	const double MAX = 1;

	int sizeA_str;
	int sizeA_stl;
	int sizeB_str;
	int sizeB_stl;
	int sizeE_str;
	int sizeE_stl;
	int sizeG_str;
	int sizeG_stl;
	int sizeC_str;
	int sizeC_stl;

	double n;
	int ti_sloj;
	int ti_vich;
	int ti_ymn;
	int num_takt_sloj;
	int num_takt_vich;
	int num_takt_ymn;

	int posled_sloj_takt;
	int posled_vich_takt;
	int posled_ymn_takt;

	double D = 0;
	int r;

	int p;
	int q;
	int m;

	double Lavr = 0;

public:
	OKMD(int p, int q, int m, int n, int ti_sloj, int ti_vich, int ti_ymn)
	{
		sizeA_str = p;
		sizeA_stl = m;

		sizeB_str = m;
		sizeB_stl = q;

		sizeE_str = 1;
		sizeE_stl = m;

		sizeG_str = p;
		sizeG_stl = q;

		sizeC_str = p;
		sizeC_stl = q;

		r = p * m + m * q + 1 * m + p * q;

		matrA = resize_matr(matrA, sizeA_str, sizeA_stl);
		matrB = resize_matr(matrB, sizeB_str, sizeB_stl);
		matrE = resize_matr(matrE, sizeE_str, sizeE_stl);
		matrG = resize_matr(matrG, sizeG_str, sizeG_stl);
		matrC = resize_matr(matrC, sizeC_str, sizeC_stl);

		matrD = resize_matr(matrD, p, q, m);
		matrF = resize_matr(matrF, p, q, m);

		matrA = generate_matr(matrA);
		matrA = rounding(matrA);
		cout << endl << endl << "MATR_A" << endl;
		print(matrA);
		matrB = generate_matr(matrB);
		matrB = rounding(matrB);
		cout << endl << endl << "MATR_B" << endl;
		print(matrB);
		matrE = generate_matr(matrE);
		matrE = rounding(matrE);
		cout << endl << endl << "MATR_E" << endl;
		print(matrE);
		matrG = generate_matr(matrG);
		matrG = rounding(matrG);
		cout << endl << endl << "MATR_G" << endl;
		print(matrG);

		this->p = p;
		this->q = q;
		this->m = m;

		this->n = n;
		num_takt_sloj = 0;
		num_takt_vich = 0;
		num_takt_ymn = 0;
		this->ti_sloj = ti_sloj;
		this->ti_vich = ti_vich;
		this->ti_ymn = ti_ymn;

	}
	vector<vector<double>> generate_matr(vector<vector<double>> matr)
	{
		for (int a = 0; a < matr.size(); a++)
		{
			for (int b = 0; b < matr[a].size(); b++)
			{
				matr[a][b] = fRand(MIN, MAX);
			}
		}
		return matr;
	}

	double fRand(double fMin, double fMax)
	{
		double f = (double)rand() / RAND_MAX;
		return fMin + f * (fMax - fMin);
	}

	vector<vector<double>> resize_matr(vector<vector<double>> matr, int str, int stl)
	{
		matr.resize(str);

		for (int b = 0; b < str; b++)
		{
			matr[b].resize(stl);
		}
		return matr;
	}

	vector<vector<vector<double>>> resize_matr(vector<vector<vector<double>>> matr, int x, int y, int z) {
		matr.resize(x);

		for (int b = 0; b < x; b++)
		{
			matr[b].resize(y);

			for (int c = 0; c < y; c++)
				matr[b][c].resize(z);
		}
		return matr;
	}

	vector<vector<double>> rounding(vector<vector<double>> mass)
	{
		for (int a = 0; a < mass.size(); a++)
		{
			for (int b = 0; b < mass[a].size(); b++)
			{
				mass[a][b] = round(mass[a][b] * 1000) / 1000;
			}
		}
		return mass;
	}

	void print(vector<vector<double>> mass)
	{
		for (int a = 0; a < mass.size(); a++)
		{
			for (int b = 0; b < mass[a].size(); b++)
			{
				cout << setw(7) << mass[a][b] << " ";
			}
			cout << endl;
		}
	}


	void start()
	{
		r = p * m * q;
		
		getMatrD();
		getMatrF();
		getMatrC();
		matrC = rounding(matrC);
		cout << endl << endl << "MATR_C" << endl;
		print(matrC);
		int T1 = posled_sloj_takt * ti_sloj + posled_vich_takt * ti_vich + posled_ymn_takt * ti_ymn;
		int Tn = num_takt_sloj * ti_sloj + num_takt_vich * ti_vich + num_takt_ymn * ti_ymn;
		Lavr /= r;
		double D = Tn / Lavr;
		cout << "TIME-" << Tn << endl;
		cout << "TIME_POSL-" << T1 << endl;
		cout << "Ky-" << T1 / double(Tn) << " e-" << T1 / double(Tn) / double(n) << " D-" << D << " r-" << r << endl;

	}

	void getMatrD()
	{
		for (int i = 0; i < p; i++)
		{
			for (int j = 0; j < q; j++)
			{
				for (int k = 0; k < m; k++)
				{
					matrD[i][j][k] = matrA[i][k] * matrB[k][j];
					posled_ymn_takt++;
					//Lavr += ti_ymn;

				}
			}
		}
		num_takt_ymn += ceil(p * q* m / n);
		Lavr += ti_ymn * p*q*m * 2;

	}

	void getMatrF()
	{
		for (int i = 0; i < p; i++)
		{
			for (int j = 0; j < q; j++)
			{
				for (int k = 0; k < m; k++)
				{
					double im = implAB(i, j, k);
					matrF[i][j][k] = im * (2 * matrE[0][k] - 1)*matrE[0][k] + implBA(i, j, k)*(1 + ((4 * im - 2)*matrE[0][k]))*(1 - matrE[0][k]);
					posled_sloj_takt += 2;
					posled_vich_takt += 1;
					if (j == 0 && i == 0) posled_vich_takt += 2;
					posled_ymn_takt += 5;
					if (j == 0 && i == 0) posled_ymn_takt += 2;
					//Lavr += 2 * ti_sloj + 3 * ti_vich + 7 * ti_ymn;

				}
			}
		}
		num_takt_sloj += ceil(p*q*m / n) * 4;

		num_takt_vich += ceil(p*m / n);
		num_takt_vich += ceil(q*m / n);
		num_takt_vich += ceil(m / n) * 2;
		num_takt_vich += ceil(p*q*m / n);

		num_takt_ymn += ceil(p*q*m / n) * 7;
		num_takt_ymn += ceil(m / n) * 2;

		Lavr += (p*m + q * m + 2 * m + p * q*m)*ti_vich + (p*q*m * 3 + p * q*m * 2)*ti_sloj + (p*q*m + m + m * 2 + p * q*m * 2 * 6)*ti_ymn;

	}

	void getMatrC()

	{
		for (int i = 0; i < p; i++)
		{
			for (int j = 0; j < q; j++)
			{
				double pr1 = prF(i, j);
				double pr2 = prD(i, j);
				matrC[i][j] = pr1 * (3 * matrG[i][j] - 2)*matrG[i][j] + (pr2 + (4 * (pr1*pr2) - 3 * pr2*matrG[i][j])*(1 - matrG[i][j]));
				posled_sloj_takt += 2;
				posled_vich_takt += 4;
				posled_ymn_takt += 8;
				//Lavr += 2 * ti_sloj + 4* ti_vich + 8 * ti_ymn;
			}
		}

		num_takt_sloj += ceil(p*q / n) * 2;

		num_takt_vich += ceil(p*q / n) * 4;
		num_takt_vich += ceil(p*q*m / n);

		num_takt_ymn += ceil(p*q / n) * 8;
		num_takt_ymn += ceil(p*q*(m - 1) / n) * 2;

		Lavr += ti_sloj * (p*q * 2 * 2) + ti_vich * (p*q*m + p * q * 3 + p * q * 2) + ti_ymn * (p*q*(m - 1) * 2 * 2 + p * q * 3 + p * q * 2 * 5);


	}

	double prF(int i, int j)
	{
		double ch = 1;
		for (int k = 0; k < m; k++)
		{
			ch *= matrF[i][j][k];
			if (k != 0)
			{
				posled_ymn_takt++;
				//num_takt_ymn++;
			}
		}
		//Lavr += (m - 1)*ti_ymn;
		return ch;

	}

	double prD(int i, int j)
	{
		double ch = 1;
		for (int k = 0; k < m; k++)
		{
			ch *= 1 - matrD[i][j][k];
			posled_vich_takt++;
			//num_takt_vich++;
			if (k != 0)
			{
				posled_ymn_takt++;
				//num_takt_ymn++;
			}

			//num_takt_vich++;
		}
		//posled_vich_takt++;
			//Lavr +=  (m-1)*ti_ymn;
		return 1 - ch;

	}

	double implAB(int i, int j, int k)
	{
		posled_sloj_takt++;
		if (i == 0)
			posled_vich_takt++;
		posled_ymn_takt++;
		return matrA[i][k] * (1 - matrB[k][j]) + 1;
	}


	double implBA(int i, int j, int k)
	{
		posled_sloj_takt++;
		if (j == 0)
			posled_vich_takt++;
		posled_ymn_takt++;
		return matrB[k][j] * (1 - matrA[i][k]) + 1;
	}
};


int main()
{
	int p, q, m, n, ti_sloj, ti_vich, ti_ymn;
	cout << "Enter p" << endl;
	cin >> p;
	cout << "Enter q" << endl;
	cin >> q;
	cout << "Enter m" << endl;
	cin >> m;
	cout << "Enter n" << endl;
	cin >> n;
	cout << "Enter ti_sloj" << endl;
	cin >> ti_sloj;
	cout << "Enter ti_vich" << endl;
	cin >> ti_vich;
	cout << "Enter ti_ymn" << endl;
	cin >> ti_ymn;

	OKMD a(p, q, m, n, ti_sloj, ti_vich, ti_ymn);
	a.start();
	system("pause");
	return 0;
}

