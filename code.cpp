#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iomanip>

#include <chrono>
#include <cstdlib>
#include <random>
#include <time.h>
using namespace std;

int deleteElement(size_t** arr, int n, int i, int x)
{
	// Search x in array
	int k;
	for (k = 2; k < n; k++)
		if (arr[i][k] == x)
			break;

	// If x found in array
	if (k < n)
	{
		// reduce size of array and move all
		// elements on space ahead
		n = n - 1;
		for (int j = k; j < n; j++)
			arr[i][j] = arr[i][j + 1];
	}
	return n;

}

struct pointcord
{
	double xcord, ycord;
};

struct point2
{
	double x, y;
	int num;
	point2* next;
};

point2* head = NULL;

void insertnode(double xcord, double ycord, int number)
{
	point2* new_node, * last;
	new_node = new point2;
	new_node->x = xcord;
	new_node->y = ycord;
	new_node->num = number;

	if (head == NULL)
	{
		head = new_node;
		new_node->next = NULL;
	}
	else
	{
		last = head;
		while (last->next != NULL)
		{
			last = last->next;
		}

		last->next = new_node;
		new_node->next = NULL;
	}
}

pointcord secondpoint(int number)
{
	pointcord cordinates;
	double a, b;
	point2* current, * previous;
	current = head;
	previous = head;
	if (current->num == number)
	{
		a = current->x;
		b = current->y;
		cordinates = { a ,b };
		return cordinates;
	}
	else
	{
		while (current->num != number)
		{
			previous = current;
			current = current->next;

		}

		a = current->x;
		b = current->y;
		cordinates = { a ,b };
		return cordinates;
	}
}

double max(double d1, double d2, double d3, double d4)
{
	if (d1 > d2 && d1 > d3 && d1 > d4) return d1;
	if (d2 > d1 && d2 > d3 && d2 > d4) return d2;
	if (d3 > d1 && d3 > d2 && d3 > d4) return d3;
	if (d4 > d1 && d4 > d2 && d4 > d3) return d4;
}

bool check(size_t n, size_t c, size_t Nvr, size_t Nhr, size_t Nv, double sr, std::vector <size_t> cells0, int* bool1, double* x_cord, double* y_cord, double r2, size_t N)
{
	double d1(0.0), d2(0.0), d3(0.0), d4(0.0), d(0.0);
	long double yc1(0.0), yc2(0.0), yc3(0.0), yc4(0.0);
	long double xc1(0.0), xc2(0.0), xc3(0.0), xc4(0.0);
	size_t* twenty1 = new size_t[20];


	twenty1[0] = n - 1 - 2 * Nv;
	twenty1[1] = n - 0 - 2 * Nv;
	twenty1[2] = n + 1 - 2 * Nv;


	twenty1[3] = n - 2 - Nv;
	twenty1[4] = n - 1 - Nv;
	twenty1[5] = n - 0 - Nv;
	twenty1[6] = n + 1 - Nv;
	twenty1[7] = n + 2 - Nv;


	twenty1[8] = n - 2;
	twenty1[9] = n - 1;
	twenty1[10] = n + 1;
	twenty1[11] = n + 2;


	twenty1[12] = n - 2 + Nv;
	twenty1[13] = n - 1 + Nv;
	twenty1[14] = n + 0 + Nv;
	twenty1[15] = n + 1 + Nv;
	twenty1[16] = n + 2 + Nv;


	twenty1[17] = n - 1 + 2 * Nv;
	twenty1[18] = n + 0 + 2 * Nv;
	twenty1[19] = n + 1 + 2 * Nv;


	xc1 = (sr * size_t(c / (Nvr)));
	xc2 = xc1;
	xc3 = xc1 + sr;
	xc4 = xc3;

	yc1 = (sr * (c % Nvr));
	yc2 = yc1 + sr;
	yc3 = yc1;
	yc4 = yc2;


	for (size_t i = 0; i < 20; i++)
	{
		if (twenty1[i] < 0 || twenty1[i] > N)
		{
			continue;
		}

		if (bool1[twenty1[i]] == 0)
		{
			continue;
		}

		if (bool1[twenty1[i]] == 2)
		{
			size_t counter(0);
			pointcord  t = secondpoint(twenty1[i]);
			d1 = (t.xcord - xc1) * (t.xcord - xc1) + (t.ycord - yc1) * (t.ycord - yc1);
			d2 = (t.xcord - xc2) * (t.xcord - xc2) + (t.ycord - yc2) * (t.ycord - yc2);
			d3 = (t.xcord - xc3) * (t.xcord - xc3) + (t.ycord - yc3) * (t.ycord - yc3);
			d4 = (t.xcord - xc4) * (t.xcord - xc4) + (t.ycord - yc4) * (t.ycord - yc4);
			d = max(d1, d2, d3, d4);

			if (d > r2)
			{
				counter++;
			}
			else
			{
				return false;
			}


			d1 = (x_cord[twenty1[i]] - xc1) * (x_cord[twenty1[i]] - xc1) + (y_cord[twenty1[i]] - yc1) * (y_cord[twenty1[i]] - yc1);
			d2 = (x_cord[twenty1[i]] - xc2) * (x_cord[twenty1[i]] - xc2) + (y_cord[twenty1[i]] - yc2) * (y_cord[twenty1[i]] - yc2);
			d3 = (x_cord[twenty1[i]] - xc3) * (x_cord[twenty1[i]] - xc3) + (y_cord[twenty1[i]] - yc3) * (y_cord[twenty1[i]] - yc3);
			d4 = (x_cord[twenty1[i]] - xc4) * (x_cord[twenty1[i]] - xc4) + (y_cord[twenty1[i]] - yc4) * (y_cord[twenty1[i]] - yc4);
			d = max(d1, d2, d3, d4);

			if (d > r2)
			{
				counter++;
			}
			else
			{
				return false;
			}
			if (counter == 2) continue;
		}
		else
		{
			// the distances between the circle center and the corneres of the cell
			d1 = (x_cord[twenty1[i]] - xc1) * (x_cord[twenty1[i]] - xc1) + (y_cord[twenty1[i]] - yc1) * (y_cord[twenty1[i]] - yc1);
			d2 = (x_cord[twenty1[i]] - xc2) * (x_cord[twenty1[i]] - xc2) + (y_cord[twenty1[i]] - yc2) * (y_cord[twenty1[i]] - yc2);
			d3 = (x_cord[twenty1[i]] - xc3) * (x_cord[twenty1[i]] - xc3) + (y_cord[twenty1[i]] - yc3) * (y_cord[twenty1[i]] - yc3);
			d4 = (x_cord[twenty1[i]] - xc4) * (x_cord[twenty1[i]] - xc4) + (y_cord[twenty1[i]] - yc4) * (y_cord[twenty1[i]] - yc4);
			d = max(d1, d2, d3, d4);

			if (d > r2)
			{
				continue;
			}
			else
			{
				return false;
			}
		}
	}

	//if (counter == 0)
	//{
	return true;
	//}
	//else return false;
}

struct cells
{
	double xs, ys, xe, ye, xd, yd;
	bool varify;
};

double delta_x(double x1, double x2, double y1, double y2, double d)
{
	double m, x;
	m = (y2 - y1) / (x2 - x1);
	x = sqrt((d * d) / (1 + (m * m)));

	if ((x2 - x1) == 0)
	{
		x = 0;
	}
	else
	{
		m = (y2 - y1) / (x2 - x1);
		x = sqrt((d * d) / (1 + (m * m)));
	}

	if ((x2 < x1))
	{
		return -x;
	}

	return x;
}

double delta_y(double x1, double x2, double y1, double y2, double deltax, double d)
{
	double y;
	double m = (y2 - y1) / (x2 - x1);

	if ((x2 - x1) == 0 && (y2 > y1))
	{
		y = d;
	}
	else if ((x2 - x1) == 0 && (y1 > y2))
	{

		y = -d;

	}
	else {

		y = m * abs(deltax);
	}

	if ((x2 < x1) && (y2 < y1))
	{
		return -y;
	}

	if ((x2 < x1) && (y2 > y1))
	{
		return -y;
	}

	return y;
}

int get_size(double* arrayx, double* arrayy, double r1, int size)
{
	int counter = 0;
	for (int i = 0; i < size; i++)
	{
		double di = sqrt((arrayx[i] - arrayx[i + 1]) * (arrayx[i] - arrayx[i + 1]) + (arrayy[i] - arrayy[i + 1]) * (arrayy[i] - arrayy[i + 1]));
		double d = sqrt((arrayx[i] - arrayx[i + 1]) * (arrayx[i] - arrayx[i + 1]) + (arrayy[i] - arrayy[i + 1]) * (arrayy[i] - arrayy[i + 1]));
		while (di > r1)
		{
			di = di / 2;
		}

		int num_segments = d / di;
		counter += (num_segments);
	}
	return (counter);
}

void cells_spilit(cells* cell, double* arrayx, double* arrayy, double r1, int size)
{
	int num_segments = 0;
	int counter = 0;

	for (int i = 0; i < size; i++)
	{
		double di = sqrt((arrayx[i] - arrayx[i + 1]) * (arrayx[i] - arrayx[i + 1]) + (arrayy[i] - arrayy[i + 1]) * (arrayy[i] - arrayy[i + 1]));
		double d = sqrt((arrayx[i] - arrayx[i + 1]) * (arrayx[i] - arrayx[i + 1]) + (arrayy[i] - arrayy[i + 1]) * (arrayy[i] - arrayy[i + 1]));

		while (di > r1)
		{
			di = di / 2;
		}

		double deltax = delta_x(arrayx[i], arrayx[i + 1], arrayy[i], arrayy[i + 1], di);
		double deltay = delta_y(arrayx[i], arrayx[i + 1], arrayy[i], arrayy[i + 1], deltax, di);

		int num_segments = d / di;

		for (size_t j = 0; j < num_segments; j++)
		{
			cell[counter].xs = arrayx[i] + j * deltax;
			cell[counter].ys = arrayy[i] + j * deltay;

			cell[counter].xe = arrayx[i] + (j + 1) * deltax;
			cell[counter].ye = arrayy[i] + (j + 1) * deltay;

			counter++;
		}
	}
}

void trimmig(size_t& no_corners, double x_cords, double y_cords, double x_cordn, double y_cordn, double* corners_x, double* corners_y, size_t* right, size_t* left, size_t* newneighbor, size_t* neighbor, size_t neighbour)
{
	bool* corners_bool = new bool[12];
	double* newcorners_x = new double[12];
	double* newcorners_y = new double[12];

	for (size_t i = 0; i < no_corners; i++)
	{
		corners_bool[i] = true;
	}

	//if (bool1[fourty_four[j]] == false) continue;
	//neighbour = fourty_four[j];
	//if (neighbour > N || neighbour < 0) continue;

	double mid_x = (x_cords + x_cordn) / 2;
	double mid_y = (y_cords + y_cordn) / 2;

	double normal_x = x_cordn - x_cords;
	double normal_y = y_cordn - y_cords;

	double V_x, V_y, validty_check;
	for (size_t k = 0; k < no_corners; k++)
	{
		V_x = corners_x[k] - mid_x;
		V_y = corners_y[k] - mid_y;
		validty_check = V_x * normal_x + V_y * normal_y;
		if (validty_check > 0)
		{
			corners_bool[k] = false;
		}
	}

	double c1_x, c1_y, c2_x, c2_y; int h = 0;
	for (size_t k = 0; k < no_corners; k++)
	{
		if ((corners_bool[k] == true) && (corners_bool[right[k]] == false))
		{
			double corner_mid_x = mid_x - corners_x[k];
			double corner_mid_y = mid_y - corners_y[k];

			double corneri_cornerj_x = corners_x[right[k]] - corners_x[k];
			double corneri_cornerj_y = corners_y[right[k]] - corners_y[k];

			double trimming_ratio = (corner_mid_x * normal_x + corner_mid_y * normal_y) / (corneri_cornerj_x * normal_x + corneri_cornerj_y * normal_y);
			c1_x = corners_x[k] + trimming_ratio * (corneri_cornerj_x);
			c1_y = corners_y[k] + trimming_ratio * (corneri_cornerj_y);
		}

		if ((corners_bool[k] == true) && (corners_bool[left[k]] == false))
		{
			h = k;
			double corner_mid_x = mid_x - corners_x[k];
			double corner_mid_y = mid_y - corners_y[k];

			double corneri_cornerj_x = corners_x[left[k]] - corners_x[k];
			double corneri_cornerj_y = corners_y[left[k]] - corners_y[k];

			double trimming_ratio = (corner_mid_x * normal_x + corner_mid_y * normal_y) / (corneri_cornerj_x * normal_x + corneri_cornerj_y * normal_y);
			c2_x = corners_x[k] + trimming_ratio * (corneri_cornerj_x);
			c2_y = corners_y[k] + trimming_ratio * (corneri_cornerj_y);
		}
	}

	int counterr(-1);
	for (size_t k = 0; k < no_corners; k++)
	{
		if ((corners_bool[k] == true) && (corners_bool[right[k]] == false))
		{
			counterr++;
			newcorners_x[counterr] = corners_x[k];
			newcorners_y[counterr] = corners_y[k];
			newneighbor[counterr] = neighbor[k];
			counterr++;
			newcorners_x[counterr] = c1_x;
			newcorners_y[counterr] = c1_y;
			newneighbor[counterr] = neighbour;
			counterr++;
			newcorners_x[counterr] = c2_x;
			newcorners_y[counterr] = c2_y;
			newneighbor[counterr] = neighbor[left[h]];
			continue;
		}
		if (corners_bool[k] == true)
		{
			counterr++;
			newcorners_x[counterr] = corners_x[k];
			newcorners_y[counterr] = corners_y[k];
			newneighbor[counterr] = neighbor[k];
		}
	}

	//delete[] corners_x;
	//delete[] corners_y;
	//delete[] right;
	//delete[] left;

	no_corners = counterr + 1;
	for (size_t k = 0; k < no_corners; k++)
	{
		corners_x[k] = newcorners_x[k];
		corners_y[k] = newcorners_y[k];
	}

	for (size_t k = 0; k < no_corners; k++)
	{
		neighbor[k] = newneighbor[k];
		neighbor[k] = newneighbor[k];
	}

	for (size_t k = 0; k < no_corners; k++)
	{
		if (k == 0)
		{
			right[k] = 1;
			left[k] = no_corners - 1;
			continue;
		}
		if (k == no_corners - 1)
		{
			right[k] = 0;
			left[k] = no_corners - 2;
			continue;
		}
		right[k] = k + 1;
		left[k] = k - 1;
	}
}

int plot_circles(std::string file_name, std::vector<double> Px, std::vector<double> Py, double r1)
{
	std::fstream file(file_name.c_str(), std::ios::out);
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale% one unit = one inch" << std::endl;

	double xmin(0), xmax(20);
	double ymin(0), ymax(10);

	double Lx(xmax - xmin);
	double Ly(ymax - ymin);

	double scale_x, scale_y, scale;
	double shift_x, shift_y;

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 1.0 - xmin * scale;
		shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
	}
	else
	{
		scale = scale_y;
		shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
		shift_y = 1.0 - ymin * scale;
	}

	file << shift_x << " " << shift_y << " translate" << std::endl;

	file << "/Courier findfont" << std::endl;
	file << "0.12 scalefont" << std::endl;
	file << "setfont" << std::endl;

	for (int i(0); i < Px.size(); i++)
	{
		file << "newpath" << std::endl;
		file << Px[i] * scale << " " << (Py[i]) * scale << " " << r1 * scale << " 0  360  arc " << std::endl;
		file << "closepath" << std::endl;
		file << "1 1 0 setrgbcolor" << std::endl;
		file << "fill" << std::endl;
		file << "0.0 setlinewidth" << std::endl; // Inactive seeds
	}

	for (int i(0); i < Px.size(); i++)
	{
		file << "newpath" << std::endl;
		file << Px[i] * scale << " " << (Py[i]) * scale << " " << r1 * scale << " 0  360  arc " << std::endl;
		file << "closepath" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl;
		file << "stroke" << std::endl;
		file << "0.0 setlinewidth" << std::endl; // Inactive seeds
	}

	for (int i(0); i < Px.size(); i++)
	{
		file << "newpath" << std::endl;
		file << Px[i] * scale << " " << (Py[i]) * scale << " " << 0.005 * scale << " 0  360  arc " << std::endl;
		file << "closepath" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl;
		file << "fill" << std::endl;
		file << "0.0 setlinewidth" << std::endl; // Inactive seeds
	}
	return 0;
}

int index(double x, double y, int r, double x_step, double y_step)
{
	int Nr, Nc, i;
	Nr = x / x_step; //cout << Nr << "  ";  // the number of the raw
	Nc = y / y_step; //cout << Nc << endl;  // the number of the column
	i = Nc + Nr * r;
	return i;
}

void dart_throw(int N, double r, cells* cell, int size, std::vector<double>& pox, std::vector<double>& poy)
{
	for (size_t i = 0; i < size; i++)
	{
		cell[i].varify = false;
		cell[i].xd = 0.0;
		cell[i].yd = 0.0;
	}

	int misses(0);

	while (misses < 10000)
	{
		int R_N;
		R_N = (rand() % (((N - 1) - 0)) + 0);

		if (cell[R_N].varify == 0)
		{
			double d = sqrt((cell[R_N].xe - cell[R_N].xs) * (cell[R_N].xe - cell[R_N].xs) + (cell[R_N].ye - cell[R_N].ys) * (cell[R_N].ye - cell[R_N].ys));

			double d_rand = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / d));

			double deltax = delta_x(cell[R_N].xs, cell[R_N].xe, cell[R_N].ys, cell[R_N].ye, d_rand);
			double deltay = delta_y(cell[R_N].xs, cell[R_N].xe, cell[R_N].ys, cell[R_N].ye, deltax, d_rand);

			double randx = cell[R_N].xs + deltax;
			double randy = cell[R_N].ys + deltay;

			int counter = 0;
			int counter1 = 0;
			int naighbour1 = R_N - 1;
			int naighbour2 = R_N + 1;

			while (true)
			{
				if (naighbour1 < 0) { break; }


				if (cell[naighbour1].varify == 1)
				{
					double d_neighbor_backword = (randx - cell[naighbour1].xd) * (randx - cell[naighbour1].xd) + (randy - cell[naighbour1].yd) * (randy - cell[naighbour1].yd);

					if (d_neighbor_backword >= r * r)
					{
						counter = 0;
						break;
					}
					else
					{
						counter = 1;
						break;
					}
				}

				double d_far = (randx - cell[naighbour1].xs) * (randx - cell[naighbour1].xs) + (randy - cell[naighbour1].ys) * (randy - cell[naighbour1].ys);

				if (d_far >= r * r)
				{
					break;
				}

				naighbour1--;
			}

			while (true)
			{
				if (naighbour2 > size - 1) { break; }

				double d_neighbor_forword = (randx - cell[naighbour2].xd) * (randx - cell[naighbour2].xd) + (randy - cell[naighbour2].yd) * (randy - cell[naighbour2].yd);

				if (cell[naighbour2].varify == 1)
				{
					if (d_neighbor_forword >= r * r)
					{
						counter1 = 0;
						break;
					}
					else
					{
						counter1 = 1;
						break;
					}
				}

				double d_far = (randx - cell[naighbour2].xe) * (randx - cell[naighbour2].xe) + (randy - cell[naighbour2].ye) * (randy - cell[naighbour2].ye);

				if (d_far >= r * r)
				{
					break;
				}

				naighbour2++;
			}

			if (counter == 0 && counter1 == 0)
			{
				cell[R_N].xd = randx;
				cell[R_N].yd = randy;
				cell[R_N].varify = true;
				pox.push_back(randx); // vector name . push_back( value )
				poy.push_back(randy);
				misses = 0;
			}

			else if (counter == 1 || counter1 == 1)
			{
				misses++, counter = 0;

			}
		}
		else
		{
			misses++;
		}
	}


	//std::cout << "pox      = " << pox.size() << std::endl;

	plot_circles("myplt.ps", pox, poy, r);
}


void dart_throw_refine(double r, cells* cell, std::vector<double>& pox, std::vector<double>& poy, std::vector<size_t>& refine0, size_t ff, int size)
{
	int misses(0);

	while (misses < 100)
	{
		if (refine0.size() == 0) break;
		int c;
		c = (rand() % (((refine0.size()) - 0)) + 0);

		int R_N = size_t(refine0[c] / ff);

		if (cell[R_N].varify == 1)
		{
			if (refine0.size() == 0) break;
			refine0.erase(refine0.begin() + c);

			//if ((refine0[c] % 2) == 0)
			//{
			//	refine0.erase(refine0.begin() + c);
			//	refine0.erase(refine0.begin() + c + 1);
			//}
			//else
			//{
			//	refine0.erase(refine0.begin() + c);
			//	refine0.erase(refine0.begin() + c - 1);
			//}
			//
			continue;
		}

		if (cell[R_N].varify == 0)
		{
			double d = sqrt((cell[R_N].xe - cell[R_N].xs) * (cell[R_N].xe - cell[R_N].xs) + (cell[R_N].ye - cell[R_N].ys) * (cell[R_N].ye - cell[R_N].ys));
			double ds = d / ff;
			double d_rand = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / ds));

			int num = (refine0[c] - R_N * ff);

			double deltax = delta_x(cell[R_N].xs, cell[R_N].xe, cell[R_N].ys, cell[R_N].ye, (ds * num));
			double deltay = delta_y(cell[R_N].xs, cell[R_N].xe, cell[R_N].ys, cell[R_N].ye, deltax, (ds * num));

			double deltax1 = delta_x(cell[R_N].xs, cell[R_N].xe, cell[R_N].ys, cell[R_N].ye, d_rand);
			double deltay1 = delta_y(cell[R_N].xs, cell[R_N].xe, cell[R_N].ys, cell[R_N].ye, deltax1, d_rand);

			double x1 = cell[R_N].xs + deltax;
			double y1 = cell[R_N].ys + deltay;

			double randx = x1 + deltax1;
			double randy = y1 + deltay1;

			int counter = 0;
			int counter1 = 0;
			int naighbour1 = R_N - 1;
			int naighbour2 = R_N + 1;

			while (true)
			{
				if (naighbour1 < 0) { break; }
				double d_neighbor_backword = (randx - cell[naighbour1].xd) * (randx - cell[naighbour1].xd) + (randy - cell[naighbour1].yd) * (randy - cell[naighbour1].yd);

				if (cell[naighbour1].varify == 1)
				{
					if (d_neighbor_backword >= r * r)
					{
						counter = 0;
						break;
					}
					else
					{
						counter = 1;
						break;
					}
				}

				double d_far = (randx - cell[naighbour1].xs) * (randx - cell[naighbour1].xs) + (randy - cell[naighbour1].ys) * (randy - cell[naighbour1].ys);

				if (d_far >= r * r)
				{
					break;
				}

				naighbour1--;
			}
			while (true)
			{
				if (naighbour2 > size - 1) { break; }
				double d_neighbor_forword = (randx - cell[naighbour2].xd) * (randx - cell[naighbour2].xd) + (randy - cell[naighbour2].yd) * (randy - cell[naighbour2].yd);

				if (cell[naighbour2].varify == 1)
				{
					if (d_neighbor_forword >= r * r)
					{
						counter1 = 0;
						break;
					}
					else
					{
						counter1 = 1;
						break;
					}
				}

				double d_far = (randx - cell[naighbour2].xe) * (randx - cell[naighbour2].xe) + (randy - cell[naighbour2].ye) * (randy - cell[naighbour2].ye);

				if (d_far >= r * r)
				{
					break;
				}

				naighbour2++;
			}

			if (counter == 0 && counter1 == 0)
			{
				cell[R_N].xd = randx;
				cell[R_N].yd = randy;
				cell[R_N].varify = true;
				pox.push_back(randx); // vector name . push_back( value )
				poy.push_back(randy);
				misses = 0;
				refine0.erase(refine0.begin() + c);

				/*if ((refine0[c] % 2) == 0)
				{
					refine0.erase(refine0.begin() + c);
					refine0.erase(refine0.begin() + c + 1);
				}
				else
				{
					refine0.erase(refine0.begin() + c);
					refine0.erase(refine0.begin() + c - 1);
				}*/
			}

			else if (counter == 1 || counter1 == 1)
			{
				misses++;
			}
		}
		else
		{
			misses++;
		}
	}


	//std::cout << "pox      = " << pox.size() << std::endl;

	//plot_circles("myplt.ps", pox, poy, r);
}


bool checkboundry(size_t N, size_t n, cells* cell, size_t g, double r, int size)
{
	if (cell[N].varify == true) return false;
	if (cell[N].varify == false)
	{
		int nr = N + 1;
		int nl = N - 1;

		if (N == 0 || N == size - 1) { goto end; }

		while (true)
		{
			if (nr == size - 1 && cell[nr].varify == false) { break; }
			if (cell[nr].varify == true) break;
			nr++;
		}
		if (nr == size - 1 && cell[nr].varify == false) { goto end; }
		while (true)
		{
			if (nl == 0 && cell[nl].varify == false) { break; }
			if (cell[nl].varify == true) break;
			nl--;
		}
		if (nl == 0 && cell[nl].varify == false) { goto end; }
		if (((cell[nl].xd - cell[nr].xd) * (cell[nl].xd - cell[nr].xd) + (cell[nl].yd - cell[nr].yd) * (cell[nl].yd - cell[nr].yd)) < (4 * r * r))
		{
			return false;
		}
	}

end:

	double d = sqrt((cell[N].xe - cell[N].xs) * (cell[N].xe - cell[N].xs) + (cell[N].ye - cell[N].ys) * (cell[N].ye - cell[N].ys));

	double ds = (d / g);

	int num = (n - N * g);

	double deltax = delta_x(cell[N].xs, cell[N].xe, cell[N].ys, cell[N].ye, (ds * num));
	double deltay = delta_y(cell[N].xs, cell[N].xe, cell[N].ys, cell[N].ye, deltax, (ds * num));


	double deltax1 = delta_x(cell[N].xs, cell[N].xe, cell[N].ys, cell[N].ye, ds);
	double deltay1 = delta_y(cell[N].xs, cell[N].xe, cell[N].ys, cell[N].ye, deltax1, ds);

	double x1 = cell[N].xs + deltax;
	double x2 = x1 + deltax1;

	double y1 = cell[N].ys + deltay;
	double y2 = y1 + deltay1;

	int counter = 0;
	int counter1 = 0;
	int naighbour1 = N - 1;
	int naighbour2 = N + 1;

	while (true)
	{
		if (naighbour1 < 0) { break; }
		if (cell[naighbour1].varify == 1)
		{
			double d_neighbor_backword = (x2 - cell[naighbour1].xd) * (x2 - cell[naighbour1].xd) + (y2 - cell[naighbour1].yd) * (y2 - cell[naighbour1].yd);

			if (d_neighbor_backword >= r * r)
			{
				counter = 0;
				break;
			}
			else
			{
				counter = 1;
				break;
			}
		}

		double d_far = (x1 - cell[naighbour1].xs) * (x1 - cell[naighbour1].xs) + (y1 - cell[naighbour1].ys) * (y1 - cell[naighbour1].ys);

		if (d_far >= r * r)
		{
			break;
		}

		naighbour1--;
	}

	while (true)
	{
		if (naighbour2 > size - 1) { break; }
		if (cell[naighbour2].varify == 1)
		{
			double d_neighbor_forword = (x1 - cell[naighbour2].xd) * (x1 - cell[naighbour2].xd) + (y1 - cell[naighbour2].yd) * (y1 - cell[naighbour2].yd);

			if (d_neighbor_forword >= r * r)
			{
				counter1 = 0;
				break;
			}
			else
			{
				counter1 = 1;
				break;
			}
		}

		double d_far = (x2 - cell[naighbour2].xe) * (x2 - cell[naighbour2].xe) + (y2 - cell[naighbour2].ye) * (y2 - cell[naighbour2].ye);

		if (d_far >= r * r)
		{
			break;
		}

		naighbour2++;
	}

	if (counter == 0 && counter1 == 0)
	{
		return true;
	}
	else if (counter == 1 || counter1 == 1)
	{
		return false;
	}
}








int main()
{
	//get the start time
	auto start = chrono::steady_clock::now();

	srand(time(NULL));
	double r = 0.1 * sqrt(2);



	double r1 = 0.7 * r;
	int size = 0;
	double* arrx = new double[5];
	double* arry = new double[5];

	arrx[0] = 0, arry[0] = 0;
	arrx[1] = 20, arry[1] = 0;
	arrx[2] = 20, arry[2] = 10;
	arrx[3] = 0, arry[3] = 10;
	arrx[4] = 0, arry[4] = 0;


	size = get_size(arrx, arry, r1, 4);

	int N1 = size;

	//std::cout << " size " << size << std::endl;

	cells* cell = new cells[N1];

	cells_spilit(cell, arrx, arry, r1, 4);

	std::vector<double>pox;
	std::vector<double>poy;

	dart_throw(N1, r1, cell, size, pox, poy);

	//for (int i = 0; i < size; i++) {

	//	std::cout << " xs     " << "  " << i << "  = " << cell[i].xs;
	//	std::cout << "           ys     " << "  " << i << "  =  " << cell[i].ys << std::endl;
	// 
	//	std::cout << " xe    " << "  " << i << "  = " << cell[i].xe;
	//	std::cout << "          ye     " << "  " << i << "  =  " << cell[i].ye << std::endl;

	//}

	//std::cout << " size " << size << std::endl;

	//std::cout << "pox      = " << pox.size() << std::endl;
	//plot_circles("myplt.ps", pox, poy, r1);


	std::vector<double>crossx;
	std::vector<double>crossy;

	pox.clear();
	poy.clear();

	for (size_t i = 0; i < size; i++)
	{
		if (cell[i].varify == false) continue;
		pox.push_back(cell[i].xd);
		poy.push_back(cell[i].yd);
	}

	for (int i(0); i < pox.size() - 1; i++)
	{
		double xm = (pox[i] + pox[i + 1]) / 2;
		double ym = (poy[i] + poy[i + 1]) / 2;
		double w = sqrt((pox[i] - pox[i + 1]) * (pox[i] - pox[i + 1]) + (poy[i] - poy[i + 1]) * (poy[i] - poy[i + 1]));
		double cos = (pox[i + 1] - pox[i]) / w;
		double sin = (poy[i + 1] - poy[i]) / w;
		if ((r1 * r1) - (.25 * w * w) < 0) continue;
		double k = sqrt((r1 * r1) - (0.25 * w * w));
		crossx.push_back(xm - (sin * k));
		crossy.push_back(ym + (cos * k));
	}

	double xm = (pox[0] + pox[pox.size() - 1]) / 2;
	double ym = (poy[0] + poy[pox.size() - 1]) / 2;
	double w = sqrt((pox[0] - pox[pox.size() - 1]) * (pox[0] - pox[pox.size() - 1]) + (poy[0] - poy[pox.size() - 1]) * (poy[0] - poy[pox.size() - 1]));
	double coss = (pox[pox.size() - 1] - pox[0]) / w;
	double sins = (poy[pox.size() - 1] - poy[0]) / w;
	double k = sqrt((r1 * r1) - (0.25 * w * w));
	crossx.push_back(xm + (sins * k));
	crossy.push_back(ym - (coss * k));

	//delete[] cell;

	pox.clear();
	poy.clear();

	/*
	double* arrxc = new double[5];
	double* arryc = new double[5];

	arrxc[0] = 9;
	arryc[0] = 4;
	arrxc[1] = 9;
	arryc[1] = 6;
	arrxc[2] = 11;
	arryc[2] = 6;
	arrxc[3] = 11;
	arryc[3] = 4;
	arrxc[4] = 9;
	arryc[4] = 4;

	int size_circle(0);

	size_circle = get_size(arrxc, arryc, r1, 4);
	*/

	double* arrxc = new double[400];
	double* arryc = new double[400];

	/*double xc = 8;
	double yc = 5;
	double rc = 2;

	for (int i = 0; i < 100; i++)
	{

		arrxc[i] = xc + rc * cos((3.6 * i) * (M_PI / 180));
		arryc[i] = yc + rc * sin((3.6 * i) * (M_PI / 180));

	}

	arrxc[100] = xc + rc * cos(0);
	arryc[100] = yc + rc * sin(0);*/


	// obtaning data from the text file
	std::fstream myfile("x.cord.txt");
	if (!myfile.is_open()) abort;

	std::string str, temp;
	int counters(-1);
	while (getline(myfile, str))
	{
		size_t start = 0;
		while (true)
		{
			size_t end = str.find_first_of('	', start);
			if (end == std::string::npos) end = str.size();
			temp = str.substr(start, end - start);
			if (isdigit(temp[0]))
			{
				counters++;
				arrxc[counters] = atof(temp.c_str());
			}
			start = end + 1;
			if (end == str.size()) break;
		}
	} // x.cord is ready



		// obtaning data from the text file
	std::fstream myfiles("y.cord.txt");
	if (!myfiles.is_open()) abort;

	std::string strr, tempp;
	int counterss(-1);
	while (getline(myfiles, strr))
	{
		size_t start = 0;
		while (true)
		{
			size_t end = strr.find_first_of('	', start);
			if (end == std::string::npos) end = strr.size();
			tempp = strr.substr(start, end - start);
			if (isdigit(tempp[0]))
			{
				counterss++;
				arryc[counterss] = atof(tempp.c_str());
			}
			start = end + 1;
			if (end == strr.size()) break;
		}
	} // y.cord is ready

	int size_circle(0);

	size_circle = get_size(arrxc, arryc, r1, counters);

	int N_circle = size_circle;

	//std::cout << " size " << size << std::endl;

	cells* cell_circle = new cells[N_circle];

	cells_spilit(cell_circle, arrxc, arryc, r1, counters);

	/*
	for (int i = 0; i < size_circle; i++) {
		std::cout << " xs     " << "  " << i << "  = " << cell_circle[i].xs;
		std::cout << "           ys     " << "  " << i << "  =  " << cell_circle[i].ys << std::endl;
		std::cout << " xe    " << "  " << i << "  = " << cell_circle[i].xe;
		std::cout << "          ye     " << "  " << i << "  =  " << cell_circle[i].ye << std::endl;
	}
	*/

	/*delete[] arrxc;
	delete[] arryc;*/

	dart_throw(N_circle, r1, cell_circle, size_circle, pox, poy);
	pox.clear();
	poy.clear();

	std::vector<size_t> refine0; std::vector<size_t> refine1;
	int r00(-1), r11(-1);
	for (size_t i = 0; i < size_circle; i++)
	{
		if (cell_circle[i].varify == false)
		{
			int nr = i + 1;
			int nl = i - 1;

			if (i == 0 || i == size_circle - 1) { refine0.push_back(i); continue; }

			while (true)
			{
				if (nr == size_circle - 1 && cell_circle[nr].varify == false) { break; }
				if (cell_circle[nr].varify == true) break;
				nr++;
			}
			if (nr == size_circle - 1 && cell_circle[nr].varify == false) { refine0.push_back(i); continue; }
			while (true)
			{
				if (nl == 0 && cell_circle[nl].varify == false) { break; }
				if (cell_circle[nl].varify == true) break;
				nl--;
			}
			if (nl == 0 && cell_circle[nl].varify == false) { refine0.push_back(i); continue; }
			if (((cell_circle[nl].xd - cell_circle[nr].xd) * (cell_circle[nl].xd - cell_circle[nr].xd) + (cell_circle[nl].yd - cell_circle[nr].yd) * (cell_circle[nl].yd - cell_circle[nr].yd)) < (4 * r1 * r1))
			{
				continue;
			}
			else
			{
				refine0.push_back(i);
			}
		}
	}
	double v(0);
	for (size_t i = 0; i < refine0.size(); i++)
	{
		if (cell_circle[refine0[i]].varify == 0) v++;
	}
	size_t ff = 2;
	size_t gg = 0;
	size_t nn = 0;
	while (true)
	{
		gg++;

		for (size_t i = 0; i < refine0.size(); i++)
		{
			refine1.push_back(refine0[i]);
		}

		refine0.clear();

		for (size_t i = 0; i < refine1.size(); i++)
		{

			if (gg == 1)
			{
				nn = refine1[i];
			}
			else
			{
				nn = size_t(refine1[i] / ff);
			}

			size_t c1, c2;

			// get children cells from current parent cells
			c1 = refine1[i] * 2;
			c2 = c1 + 1;


			if (checkboundry(nn, c1, cell_circle, ff, r1, N_circle)) refine0.push_back(c1);
			if (checkboundry(nn, c2, cell_circle, ff, r1, N_circle)) refine0.push_back(c2);
		}

		std::cout << refine0.size() << std::endl;
		if (refine0.size() == 0) { break; }
		/*	if ( gg == 15) { break; }*/

		dart_throw_refine(r1, cell_circle, pox, poy, refine0, ff, N_circle);

		ff = 2 * ff;
		refine1.clear();
	}

	pox.clear();
	poy.clear();

	std::vector<double>crossx1;
	std::vector<double>crossy1;

	for (size_t i = 0; i < size_circle; i++)
	{
		if (cell_circle[i].varify == false) continue;
		pox.push_back(cell_circle[i].xd);
		poy.push_back(cell_circle[i].yd);
	}

	for (int i(0); i < pox.size() - 1; i++)
	{
		double xm = (pox[i] + pox[i + 1]) / 2;
		double ym = (poy[i] + poy[i + 1]) / 2;
		double w = sqrt((pox[i] - pox[i + 1]) * (pox[i] - pox[i + 1]) + (poy[i] - poy[i + 1]) * (poy[i] - poy[i + 1]));
		double cos = (pox[i + 1] - pox[i]) / w;
		double sin = (poy[i + 1] - poy[i]) / w;
		if ((r1 * r1) - (.25 * w * w) < 0) { continue; }
		double k = sqrt((r1 * r1) - (0.25 * w * w));
		crossx1.push_back(xm + (sin * k));
		crossy1.push_back(ym - (cos * k));
		crossx1.push_back(xm - (sin * k));
		crossy1.push_back(ym + (cos * k));
	}

	xm = (pox[0] + pox[pox.size() - 1]) / 2;
	ym = (poy[0] + poy[pox.size() - 1]) / 2;
	w = sqrt((pox[0] - pox[pox.size() - 1]) * (pox[0] - pox[pox.size() - 1]) + (poy[0] - poy[pox.size() - 1]) * (poy[0] - poy[pox.size() - 1]));
	double cos = (pox[pox.size() - 1] - pox[0]) / w;
	double sin = (poy[pox.size() - 1] - poy[0]) / w;
	k = sqrt((r1 * r1) - (0.25 * w * w));
	crossx1.push_back(xm + (sin * k));
	crossy1.push_back(ym - (cos * k));
	crossx1.push_back(xm - (sin * k));
	crossy1.push_back(ym + (cos * k));

	delete[] cell_circle;






	long double r2 = r * r;
	long double s = r / sqrt(2);
	size_t Nv = size_t(10 / s);
	size_t Nh = size_t(20 / s);
	int N = Nv * Nh;

	int* bool1 = new int[N];
	double* x_cord = new double[N];
	double* y_cord = new double[N];

	for (size_t i = 0; i < N; i++)
	{

		bool1[i] = 0;
		x_cord[i] = 0.0;  //check
		y_cord[i] = 0.0;  //check
	}

	std::vector <size_t> covered;

	for (int i(0); i < crossx1.size(); i++)
	{
		int u = index(crossx1[i], crossy1[i], Nv, s, s);
		if (bool1[u] == 1)
		{
			bool1[u] = 2;
			insertnode(crossx1[i], crossy1[i], u);
		}
		else
		{

			bool1[u] = 1;
			x_cord[u] = crossx1[i];
			y_cord[u] = crossy1[i];
			covered.push_back(u);
		}
	}
	int a = covered.size();
	for (int i(0); i < crossx.size(); i++)
	{
		int u = index(crossx[i], crossy[i], Nv, s, s);
		if (bool1[u] == 1)
		{
			bool1[u] = 2;
			insertnode(crossx[i], crossy[i], u);
		}
		else
		{

			bool1[u] = 1;
			x_cord[u] = crossx[i];
			y_cord[u] = crossy[i];
			covered.push_back(u);
		}
	}

	
	/*crossx1.clear();
	crossx.clear();
	crossy1.clear();
	crossy.clear();*/
	pox.clear();
	poy.clear();

	std::vector <size_t> cells0, cells1;
	int bug = 0;

	long double x(0.0), y(0.0), X(0.0), Y(0.0);

	int* twenty = new int[20];
	int n = 0;
	size_t miss = 0;
	while (true)
	{
		static default_random_engine e;
		static uniform_int_distribution<unsigned> u(0, N - 1);
		n = u(e);

		if (bool1[n] != 0) { continue; }
		long double x = ((long double)rand() / RAND_MAX) * s;
		long double y = ((long double)rand() / RAND_MAX) * s;

		twenty[0] = n - 1 - 2 * Nv;
		twenty[1] = n - 0 - 2 * Nv;
		twenty[2] = n + 1 - 2 * Nv;


		twenty[3] = n - 2 - Nv;
		twenty[4] = n - 1 - Nv;
		twenty[5] = n - 0 - Nv;
		twenty[6] = n + 1 - Nv;
		twenty[7] = n + 2 - Nv;


		twenty[8] = n - 2;
		twenty[9] = n - 1;
		twenty[10] = n + 1;
		twenty[11] = n + 2;


		twenty[12] = n - 2 + Nv;
		twenty[13] = n - 1 + Nv;
		twenty[14] = n + 0 + Nv;
		twenty[15] = n + 1 + Nv;
		twenty[16] = n + 2 + Nv;


		twenty[17] = n - 1 + 2 * Nv;
		twenty[18] = n + 0 + 2 * Nv;
		twenty[19] = n + 1 + 2 * Nv;


		X = (s * size_t(n / Nv)) + x;
		Y = (s * (n % Nv)) + y;

		size_t check = 0;
		for (size_t i = 0; i < 20; i++)
		{
			if (twenty[i] < 0 || twenty[i] > N)
			{
				continue;
			}

			if (bool1[twenty[i]] == 0)
			{
				continue;
			}

			if (bool1[twenty[i]] == 2)
			{
				pointcord  l = secondpoint(twenty[i]);
				if ((((x_cord[twenty[i]] - X) * (x_cord[twenty[i]] - X)) + ((y_cord[twenty[i]] - Y) * (y_cord[twenty[i]] - Y))) < r2)
				{
					check++;
					miss++;
					break;
				}

				if ((((l.xcord - X) * (l.xcord - X)) + ((l.ycord - Y) * (l.ycord - Y))) < r2)
				{
					check++;
					miss++;
					break;
				}
			}
			else
			{
				if ((((x_cord[twenty[i]] - X) * (x_cord[twenty[i]] - X)) + ((y_cord[twenty[i]] - Y) * (y_cord[twenty[i]] - Y))) < r2)
				{
					check++;
					miss++;
					break;
				}
			}
		}

		if (check == 0) // point is valid
		{
			bool1[n] = 1;
			x_cord[n] = X;
			y_cord[n] = Y;
			covered.push_back(n);
			miss = 0;
		}

		if (miss > 100) break;
	}

	for (size_t i = 0; i < N; i++)
	{
		if (bool1[i] == 0) cells0.push_back(i);
	}

	size_t Nr = 4 * N, Nvr = 2 * Nv, Nhr = 2 * Nh, f = 2;
	long double sr = 0.5 * s;

	size_t g = 0;
	while (true)
	{
		g++;
		for (size_t i = 0; i < cells0.size(); i++)
		{
			cells1.push_back(cells0[i]);
		}

		cells0.clear();

		for (size_t i = 0; i < cells1.size(); i++)    // full the cells vector
		{
			if (g == 1)
			{
				n = cells1[i];
			}
			else
			{
				n = size_t(size_t(cells1[i] / (Nvr / 2)) / (f / 2)) * Nv + size_t((cells1[i] % (Nvr / 2)) / (f / 2));
			}

			size_t c1, c2, c3, c4;

			// get children cells from current parent cells
			c1 = (size_t(cells1[i] / (Nvr / 2)) * Nvr * 2) + (cells1[i] % (Nvr / 2)) * 2;
			c2 = c1 + 1;
			c3 = c1 + Nvr;
			c4 = c3 + 1;

			if (check(n, c1, Nvr, Nhr, Nv, sr, cells0, bool1, x_cord, y_cord, r2, N)) cells0.push_back(c1);
			if (check(n, c2, Nvr, Nhr, Nv, sr, cells0, bool1, x_cord, y_cord, r2, N)) cells0.push_back(c2);
			if (check(n, c3, Nvr, Nhr, Nv, sr, cells0, bool1, x_cord, y_cord, r2, N)) cells0.push_back(c3);
			if (check(n, c4, Nvr, Nhr, Nv, sr, cells0, bool1, x_cord, y_cord, r2, N)) cells0.push_back(c4);
		}

		std::cout << cells0.size() << std::endl;
		if (cells0.size() == 0) { break; }

		int* twenty2 = new int[20];
		size_t miss1 = 0;

		// dart throwing on cells vector
		while (true)
		{
			int m, n;

			m = ((double)rand() / RAND_MAX) * (cells0.size() - 1);

			/*static default_random_engine e;
			static uniform_int_distribution <unsigned> u(1, (cells.size() - 1));
			m = u(e);*/

			n = size_t(size_t(cells0[m] / (Nvr)) / (f)) * Nv + size_t((cells0[m] % (Nvr)) / (f)); // no of father cell

			if (bool1[n] != 0)
			{
				cells0.erase(cells0.begin() + m);
				if (cells0.size() == 0) break;
				continue;
			}

			long double x1 = ((double)rand() / RAND_MAX) * sr;
			long double y1 = ((double)rand() / RAND_MAX) * sr;

			twenty2[0] = n - 1 - 2 * Nv;
			twenty2[1] = n - 0 - 2 * Nv;
			twenty2[2] = n + 1 - 2 * Nv;


			twenty2[3] = n - 2 - Nv;
			twenty2[4] = n - 1 - Nv;
			twenty2[5] = n - 0 - Nv;
			twenty2[6] = n + 1 - Nv;
			twenty2[7] = n + 2 - Nv;


			twenty2[8] = n - 2;
			twenty2[9] = n - 1;
			twenty2[10] = n + 1;
			twenty2[11] = n + 2;


			twenty2[12] = n - 2 + Nv;
			twenty2[13] = n - 1 + Nv;
			twenty2[14] = n + 0 + Nv;
			twenty2[15] = n + 1 + Nv;
			twenty2[16] = n + 2 + Nv;


			twenty2[17] = n - 1 + 2 * Nv;
			twenty2[18] = n + 0 + 2 * Nv;
			twenty2[19] = n + 1 + 2 * Nv;


			double X1 = (sr * size_t(cells0[m] / Nvr)) + x1;
			double Y1 = (sr * (cells0[m] % Nvr)) + y1;

			size_t check1 = 0;
			for (size_t i = 0; i < 20; i++)
			{
				if (twenty2[i] < 0 || twenty2[i] > N)
				{
					continue;
				}

				if (bool1[twenty2[i]] == 0)
				{
					continue;
				}

				if (bool1[twenty[i]] == 2)
				{
					pointcord  l = secondpoint(twenty[i]);
					if ((((x_cord[twenty[i]] - X1) * (x_cord[twenty[i]] - X1)) + ((y_cord[twenty[i]] - Y1) * (y_cord[twenty[i]] - Y1))) < r2)
					{
						check1++;
						miss1++;
						break;
					}

					if ((((l.xcord - X1) * (l.xcord - X1)) + ((l.ycord - Y1) * (l.ycord - Y1))) < r2)
					{
						check1++;
						miss1++;
						break;
					}
				}
				else
				{
					if ((((x_cord[twenty2[i]] - X1) * (x_cord[twenty2[i]] - X1)) + ((y_cord[twenty2[i]] - Y1) * (y_cord[twenty2[i]] - Y1))) < r2)
					{
						check1++;
						miss1++;
						break;
					}
				}
			}

			if (check1 == 0) // point is valid
			{
				bool1[n] = 1;
				x_cord[n] = X1;
				y_cord[n] = Y1;
				covered.push_back(n);
				miss1 = 0;
				cells0.erase(cells0.begin() + m);
			}

			if (miss1 > 100) break;
		}

		Nvr = 2 * Nvr; Nhr = 2 * Nhr; sr = 0.5 * sr; Nr = 4 * Nr; f = 2 * f;
		cells1.clear();
	}

	//std::cout << covered.size() << std::endl;

	std::string file_name = "myplot.ps";
	std::fstream file(file_name.c_str(), std::ios::out);
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale% one unit = one inch" << std::endl;
	double xmin(0), xmax(20);
	double ymin(0), ymax(10);
	double Lx(xmax - xmin);
	double Ly(ymax - ymin);
	double scale_x, scale_y, scale;
	double shift_x, shift_y;
	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 1.0 - xmin * scale;
		shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
	}
	else
	{
		scale = scale_y;
		shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
		shift_y = 1.0 - ymin * scale;
	}
	file << shift_x << " " << shift_y << " translate" << std::endl;
	file << "/Courier findfont" << std::endl;
	file << "0.12 scalefont" << std::endl;
	file << "setfont" << std::endl;






	//string  fileName;
	//cout << "ENTER FILE NAME (**.txt) ";
	//cin >> fileName;
	//ofstream txtOut;
	//txtOut.open(fileName);
	size_t** Graph = new size_t * [covered.size()];
	double** cornersx = new double* [covered.size()];
	double** cornersy = new double* [covered.size()];


	int graphpointer = 0;

	double mid_x, mid_y, normal_x, normal_y, V_x, V_y;
	double validty_check, corner_mid_x, corner_mid_y, trimming_ratio, c1_x, c1_y, c2_x, c2_y, corneri_cornerj_x, corneri_cornerj_y;
	for (size_t i = 0; i < covered.size(); i++)
	{
		size_t* fourty_four = new size_t[44];
		double* corners_x = new double[12];
		double* corners_y = new double[12];
		double* newcorners_x = new double[12];
		double* newcorners_y = new double[12];
		size_t* neighbor = new size_t[12];
		size_t* newneighbor = new size_t[12];

		size_t* right = new size_t[12];
		size_t* left = new size_t[12];
		bool* corners_bool = new bool[12];
		size_t no_corners = 4;
		size_t neighbour, seed;

		corners_x[0] = 0;
		corners_x[1] = 20;
		corners_x[2] = 20;
		corners_x[3] = 0;

		corners_y[0] = 0;
		corners_y[1] = 0;
		corners_y[2] = 10;
		corners_y[3] = 10;

		right[0] = 1;
		right[1] = 2;
		right[2] = 3;
		right[3] = 0;

		left[0] = 3;
		left[1] = 0;
		left[2] = 1;
		left[3] = 2;

		neighbor[0] = 0;
		neighbor[1] = 0;
		neighbor[2] = 0;
		neighbor[3] = 0;


		seed = covered[i];
		fourty_four[0] = seed - 2 - 3 * Nv;
		fourty_four[1] = seed - 1 - 3 * Nv;
		fourty_four[2] = seed - 0 - 3 * Nv;
		fourty_four[3] = seed + 1 - 3 * Nv;
		fourty_four[4] = seed + 2 - 3 * Nv;


		fourty_four[5] = seed - 3 - 2 * Nv;
		fourty_four[6] = seed - 2 - 2 * Nv;
		fourty_four[7] = seed - 1 - 2 * Nv;
		fourty_four[8] = seed - 0 - 2 * Nv;
		fourty_four[9] = seed + 1 - 2 * Nv;
		fourty_four[10] = seed + 2 - 2 * Nv;
		fourty_four[11] = seed + 3 - 2 * Nv;


		fourty_four[12] = seed - 3 - Nv;
		fourty_four[13] = seed - 2 - Nv;
		fourty_four[14] = seed - 1 - Nv;
		fourty_four[15] = seed - 0 - Nv;
		fourty_four[16] = seed + 1 - Nv;
		fourty_four[17] = seed + 2 - Nv;
		fourty_four[18] = seed + 3 - Nv;


		fourty_four[19] = seed - 3;
		fourty_four[20] = seed - 2;
		fourty_four[21] = seed - 1;
		fourty_four[22] = seed + 1;
		fourty_four[23] = seed + 2;
		fourty_four[24] = seed + 3;


		fourty_four[25] = seed - 3 + Nv;
		fourty_four[26] = seed - 2 + Nv;
		fourty_four[27] = seed - 1 + Nv;
		fourty_four[28] = seed + 0 + Nv;
		fourty_four[29] = seed + 1 + Nv;
		fourty_four[30] = seed + 2 + Nv;
		fourty_four[31] = seed + 3 + Nv;


		fourty_four[32] = seed - 3 + 2 * Nv;
		fourty_four[33] = seed - 2 + 2 * Nv;
		fourty_four[34] = seed - 1 + 2 * Nv;
		fourty_four[35] = seed + 0 + 2 * Nv;
		fourty_four[36] = seed + 1 + 2 * Nv;
		fourty_four[37] = seed + 2 + 2 * Nv;
		fourty_four[38] = seed + 3 + 2 * Nv;


		fourty_four[39] = seed - 2 + 3 * Nv;
		fourty_four[40] = seed - 1 + 3 * Nv;
		fourty_four[41] = seed - 0 + 3 * Nv;
		fourty_four[42] = seed + 1 + 3 * Nv;
		fourty_four[43] = seed + 2 + 3 * Nv;

		if (bool1[seed] == 2)
		{
			pointcord  c = secondpoint(seed);
			neighbour = seed;
			trimmig(no_corners, x_cord[seed], y_cord[seed], c.xcord, c.ycord, corners_x, corners_y, right, left, newneighbor, neighbor, neighbour);

			for (size_t j = 0; j < 44; j++)
			{
				if (bool1[fourty_four[j]] == 0) continue;

				neighbour = fourty_four[j];
				if (neighbour > N || neighbour < 0) continue;

				if (bool1[neighbour] == 2)
				{
					pointcord  o = secondpoint(neighbour);
					trimmig(no_corners, x_cord[seed], y_cord[seed], o.xcord, o.ycord, corners_x, corners_y, right, left, newneighbor, neighbor, 0);
					trimmig(no_corners, x_cord[seed], y_cord[seed], x_cord[neighbour], y_cord[neighbour], corners_x, corners_y, right, left, newneighbor, neighbor, 0);
				}
				else
				{
					trimmig(no_corners, x_cord[seed], y_cord[seed], x_cord[neighbour], y_cord[neighbour], corners_x, corners_y, right, left, newneighbor, neighbor, neighbour);
				}
			}
		}
		else
		{
			for (size_t j = 0; j < 44; j++)
			{
				if (bool1[fourty_four[j]] == 0) continue;

				for (size_t i = 0; i < no_corners; i++)
				{
					corners_bool[i] = true;
				}

				neighbour = fourty_four[j];
				if (neighbour > N || neighbour < 0) continue;

				if (bool1[neighbour] == 2)
				{
					pointcord  f = secondpoint(neighbour);
					trimmig(no_corners, x_cord[seed], y_cord[seed], f.xcord, f.ycord, corners_x, corners_y, right, left, newneighbor, neighbor, 0);
					trimmig(no_corners, x_cord[seed], y_cord[seed], x_cord[neighbour], y_cord[neighbour], corners_x, corners_y, right, left, newneighbor, neighbor, 0);
				}
				else
				{
					trimmig(no_corners, x_cord[seed], y_cord[seed], x_cord[neighbour], y_cord[neighbour], corners_x, corners_y, right, left, newneighbor, neighbor, neighbour);
				}
			}
		}

		Graph[graphpointer] = new size_t[no_corners + 2];
		Graph[graphpointer][0] = seed;
		Graph[graphpointer][1] = no_corners;

		for (size_t k = 0; k < no_corners; k++)
		{
			Graph[graphpointer][k + 2] = neighbor[k];
		}

		cornersx[graphpointer] = new double[no_corners];
		for (size_t k = 0; k < no_corners; k++)
		{
			cornersx[graphpointer][k] = corners_x[k];
		}

		cornersy[graphpointer] = new double[no_corners];
		for (size_t k = 0; k < no_corners; k++)
		{
			cornersy[graphpointer][k] = corners_y[k];
		}



		graphpointer++;

	}

	for (int i(0); i < a; i++)
	{
		for (int k(0); k < a; k++)
		{
			for (int o(2); o < (Graph[k][1]+2); o++)
			{
				/*if (bool1[Graph[k][o]] == 2) Graph[k][o] = 0;*/
				if (covered[i] == Graph[k][o])
				{
					Graph[k][o] = 0;
				}
			}
		}
	}

	/*for (int i = 0; i < covered.size(); i++)
	{
		for (int j = 0; j < (Graph[i][1] + 2 ); j++)
		{

			cout << Graph[i][j] << "\t";

		}
		cout << endl;
	}*/

	size_t* sc = new size_t[N];

	for (size_t i = 0; i < covered.size(); i++)
	{
		sc[Graph[i][0]] = 0;
	}

	for (int i(0); i < a; i++)
	{
		sc[Graph[i][0]] = 512;
	}

	size_t cc(1);
	while (true)
	{
		size_t fc;
		size_t count(0);

		for (size_t i = 0; i < covered.size(); i++)
		{
			if (sc[Graph[i][0]] == 0) { fc = Graph[i][0]; count++; break; }
		}

		/*if (cc == 3)  break;*/
		if(count == 0) break;
		sc[fc] = cc;

		while (true)
		{
			bool done = true;
			for (size_t i = 0; i < covered.size(); i++)
			{
				/*if (bool1[Graph[i][0]] == 2) { sc[Graph[i][0]] = 50;  continue; }*/
				if (sc[Graph[i][0]] == cc)
				{
					for (int k(2); k < (Graph[i][1] + 2); k++)
					{
						if (Graph[i][k] == 0) { sc[Graph[i][k]] = 500;  continue; }
						/*if (bool1[Graph[i][k]] == 2) { sc[Graph[i][k]] = 5000;  continue; }*/
						if (sc[Graph[i][k]] == 0)
						{
							sc[Graph[i][k]] = cc;
							done = false;


						}
					}
				}
			}
			if (done) break;
		}
		cc++;
	}


	for (size_t i = 0; i < covered.size(); i++)
	{
		if (sc[Graph[i][0]] == 1)
		{
			for (size_t k = 0; k < (Graph[i][1]-1); k++)
			{
				file << "newpath" << std::endl;
				file << cornersx[i][k] * scale << " " << (cornersy[i][k]) * scale << " ";
				file << "moveto" << std::endl;
				file << cornersx[i][k+1] * scale << " " << (cornersy[i][k+1]) * scale << " ";
				file << "lineto" << std::endl;
				file << ".003 setlinewidth" << std::endl; // Inactive seeds
				file << "1 0 0 setrgbcolor" << std::endl;
				file << "stroke" << std::endl;
			}
			file << "newpath" << std::endl;
			file << cornersx[i][0] * scale << " " << (cornersy[i][0]) * scale << " ";
			file << "moveto" << std::endl;
			file << cornersx[i][Graph[i][1]-1] * scale << " " << (cornersy[i][Graph[i][1]-1]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".003 setlinewidth" << std::endl; // Inactive seeds
			file << "1 0 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;
		}

		if (sc[Graph[i][0]] == 2)
		{
			for (size_t k = 0; k < (Graph[i][1] - 1); k++)
			{
				file << "newpath" << std::endl;
				file << cornersx[i][k] * scale << " " << (cornersy[i][k]) * scale << " ";
				file << "moveto" << std::endl;
				file << cornersx[i][k+1] * scale << " " << (cornersy[i][k+1]) * scale << " ";
				file << "lineto" << std::endl;
				file << ".003 setlinewidth" << std::endl; // Inactive seeds
				file << "0 1 0 setrgbcolor" << std::endl;
				file << "stroke" << std::endl;
			}
			file << "newpath" << std::endl;
			file << cornersx[i][0] * scale << " " << (cornersy[i][0]) * scale << " ";
			file << "moveto" << std::endl;
			file << cornersx[i][Graph[i][1]-1] * scale << " " << (cornersy[i][Graph[i][1]-1]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".003 setlinewidth" << std::endl; // Inactive seeds
			file << "0 1 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;
		}

		if (sc[Graph[i][0]] == 512)
		{
			for (size_t k = 0; k < (Graph[i][1] - 1); k++)
			{
				file << "newpath" << std::endl;
				file << cornersx[i][k] * scale << " " << (cornersy[i][k]) * scale << " ";
				file << "moveto" << std::endl;
				file << cornersx[i][k + 1] * scale << " " << (cornersy[i][k + 1]) * scale << " ";
				file << "lineto" << std::endl;
				file << ".003 setlinewidth" << std::endl; // Inactive seeds
				file << "0 0 1 setrgbcolor" << std::endl;
				file << "stroke" << std::endl;
			}
			file << "newpath" << std::endl;
			file << cornersx[i][0] * scale << " " << (cornersy[i][0]) * scale << " ";
			file << "moveto" << std::endl;
			file << cornersx[i][Graph[i][1] - 1] * scale << " " << (cornersy[i][Graph[i][1] - 1]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".003 setlinewidth" << std::endl; // Inactive seeds
			file << "0 0 1 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;
		}

		if (sc[Graph[i][0]] == 3)
		{
			for (size_t k = 0; k < (Graph[i][1] - 1); k++)
			{
				file << "newpath" << std::endl;
				file << cornersx[i][k] * scale << " " << (cornersy[i][k]) * scale << " ";
				file << "moveto" << std::endl;
				file << cornersx[i][k + 1] * scale << " " << (cornersy[i][k + 1]) * scale << " ";
				file << "lineto" << std::endl;
				file << ".003 setlinewidth" << std::endl; // Inactive seeds
				file << "0 5 0 setrgbcolor" << std::endl;
				file << "stroke" << std::endl;
			}
			file << "newpath" << std::endl;
			file << cornersx[i][0] * scale << " " << (cornersy[i][0]) * scale << " ";
			file << "moveto" << std::endl;
			file << cornersx[i][Graph[i][1] - 1] * scale << " " << (cornersy[i][Graph[i][1] - 1]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".003 setlinewidth" << std::endl; // Inactive seeds
			file << "0 5 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;
		}

		if (sc[Graph[i][0]] == 4)
		{
			for (size_t k = 0; k < (Graph[i][1] - 1); k++)
			{
				file << "newpath" << std::endl;
				file << cornersx[i][k] * scale << " " << (cornersy[i][k]) * scale << " ";
				file << "moveto" << std::endl;
				file << cornersx[i][k + 1] * scale << " " << (cornersy[i][k + 1]) * scale << " ";
				file << "lineto" << std::endl;
				file << ".003 setlinewidth" << std::endl; // Inactive seeds
				file << "0 8 0 setrgbcolor" << std::endl;
				file << "stroke" << std::endl;
			}
			file << "newpath" << std::endl;
			file << cornersx[i][0] * scale << " " << (cornersy[i][0]) * scale << " ";
			file << "moveto" << std::endl;
			file << cornersx[i][Graph[i][1] - 1] * scale << " " << (cornersy[i][Graph[i][1] - 1]) * scale << " ";
			file << "lineto" << std::endl;
			file << ".003 setlinewidth" << std::endl; // Inactive seeds
			file << "0 8 0 setrgbcolor" << std::endl;
			file << "stroke" << std::endl;
		}
	}


	// get the end time 
	auto end = chrono::steady_clock::now();

	//find the difference
	double elapsed_time_ns = double(std::chrono::duration_cast <std::chrono::nanoseconds> (end - start).count());

	// output
	std::cout << "elapsed time(s): " << elapsed_time_ns / 1e9 << endl;

}















