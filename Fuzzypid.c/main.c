#include <iostream>
#include "FuzzyPID.h"
int main()
{
    FuzzyPID myfuzzypid;
    float Target = 600;
    float actual = 0;
    float e_max =1000;
    float e_min = -1000;
    float ec_max = 800;
    float ec_min = -800;
    float kp_max =100;
    float kp_min = -100;
    float ki_max = 0.1;
    float ki_min = -0.1;
    float kd_max = 0.01;
    float kd_min = -0.01;
    float erro;
    float erro_c;
    float erro_pre = 0;
    float erro_ppre = 0;
    erro =Target - actual;
    erro_c = erro - erro_pre;
    for (int i = 0; i < 100; i++)
    {
        float u;
        u = myfuzzypid.FuzzyPIDcontroller(e_max, e_min, ec_max, ec_min, kp_max, kp_min, erro, erro_c,ki_max,ki_min,kd_max,kd_min,erro_pre,erro_ppre);
        actual +=u;
        erro_ppre = erro_pre;
        erro_pre = erro;
        erro = Target - actual;
        erro_c= erro - erro_pre;
        std::cout << "i:" << i << "\t" << "Target:" << Target << "\t" << "Actual:" << actual  << std::endl;
    }
}
