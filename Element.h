#include "Node.h"
class Element {
	Node node1;
	Node node2;
	Node node3;
	Node node4;
	int *ID;
	float **jakobian;
	float *detJ;
	float **jakobian1;
	float **matrixH;
	float **matrixC;
	float **matrixHBC;
	float *P;
	bool flagal = 0;
	bool flagap = 0;
	bool flagag = 0;
	bool flagad = 0;
	float initial_temperature;
	float simulation_time;
	float simulation_step_time;
	float To;
	float alfa;
	float c;
	float k; 
	float ro;
	float H;
	float B;
	int N_H;
	int N_B;
public:
	void fun();
	void fun1(float***&, float**&, int, int);
	float ***fun2(float***);
	void det(float***&);
	float mnoz();
	void matrix();
	void matrixc();
	void matrixhbc();
	void p();
	int get_id(int);
	int get_p(int);
	float get_matrixh(int, int);
	float get_matrixc(int, int);
	float get_matrixhbc(int, int);
	float get_nh();
	float get_nb();
	void set_detJ(float*);
	void set_Jakobian(float**);
	void set_Jakobian1(float**);
	void set_matrixH(float**);
	void set_matrixC(float**);
	void set_matrixHBC(float**);
	void set_flagal(bool);
	void set_flagap(bool);
	void set_flagag(bool);
	void set_flagad(bool);
	Element();
	void set(int, Node*, float, float, float, float, float, float, float, int, int, float, float, float);
	~Element();
};