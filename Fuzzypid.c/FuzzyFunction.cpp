/**
 * @file      FuzzyFunction.cpp/.c
 * @brief     Fuzzy PID controller implementation for embedded systems
 * @author    Xie Bingwei
 * @date      2025.12.30
 * @version   v1.1.0
 * @license   MIT License
 */
#include "FuzzyPID.h"
FuzzyPID::FuzzyPID()
{
    kp = 0;
    ki = 0;
    kd = 0;
    fuzzy_output = 0;
    qdetail_kp = 0;
    qdetail_ki = 0;
    qdetail_kd = 0;
    qfuzzy_output = 0;
    errosum = 0;
}

FuzzyPID::~FuzzyPID() {}

/**
 * @brief     This module implements the membership degree calculation function of the fuzzy PID controller.
 *            Through the triangular membership function, the precise error and error change rate are mapped
 *            to the membership degree of the fuzzy set.
 * @param[in] erro: Current system error
 * @param[in] erro_c: Error change rate (de/dt ≈ error - error_last)
 * @return    void
 */
void FuzzyPID::Get_grad_membership(float erro,float erro_c)
{
    if (erro > e_membership_values[0] && erro < e_membership_values[6])
    {
        for (int i = 0; i < num_area - 2; i++) // 因为需要检测6个区间->[-3,-2]...[2,3],故仅需要6次循环
        {
            if (erro >= e_membership_values[i] && erro <= e_membership_values[i + 1])
            {
                e_gradmembership[0] = (e_membership_values[i + 1] - erro) / (e_membership_values[i + 1] - e_membership_values[i]); // 对左边模糊集(PS)的隶属度->距右边界距离/区间长度
                e_gradmembership[1] = (erro - e_membership_values[i]) / (e_membership_values[i + 1] - e_membership_values[i]); // 对右边模糊集(PM)的隶属度->距左边界的距离/区间长度
                e_grad_index[0] = i;
                e_grad_index[1] = i + 1;  // 为对应论域区间的两个索引
                break;
            }
        }
    }
    else
    {
        if (erro <= e_membership_values[0]) // 如果误差比最小的隶属值(即-3)还小，则其完全属于NB
        {
            e_gradmembership[0] = 1; // 超出最小隶属值部分默认为0,则此时err直接等于-3,则依据算法(如上)左隶属度为1,右为0
            e_gradmembership[1] = 0;
            e_grad_index[0] = 0; // 0为NB的索引
            e_grad_index[1] = -1; // -1作为无效索引
        }
        else if (erro >= e_membership_values[6]) // 如果差值比最大的隶属值(即3)还大，则完全属于PB,err值直接等于3
        {
            e_gradmembership[0] = 1;
            e_gradmembership[1] = 0;
            e_grad_index[0] = 6; // 6为PB的索引
            e_grad_index[1] = -1; // 同上
        }
    }
    // 以下为de/dt的情况
    if (erro_c > ec_membership_values[0] && erro_c < ec_membership_values[6]) //这一部分是对de/dt隶属值的计算，方法同上
    {
        for (int i = 0; i < num_area - 2; i++)
        {
            if (erro_c >= ec_membership_values[i] && erro_c <=  ec_membership_values[i + 1])
            {
                ec_gradmembership[0] = (ec_membership_values[i + 1] - erro_c) / (ec_membership_values[i + 1] - ec_membership_values[i]);
                ec_gradmembership[1] = (erro_c - ec_membership_values[i]) / (ec_membership_values[i + 1] - ec_membership_values[i]);
                ec_grad_index[0] = i;
                ec_grad_index[1] = i + 1;
                break;
            }
        }
    }
    else
    {
        if (erro_c <= ec_membership_values[0])
        {
            ec_gradmembership[0] = 1;
            ec_gradmembership[1] = 0;
            ec_grad_index[0] = 0;
            ec_grad_index[1] = -1;
        }
        else if (erro_c >= ec_membership_values[6])
        {
            ec_gradmembership[0] = 1;
            ec_gradmembership[1] = 0;
            ec_grad_index[0] = 6;
            ec_grad_index[1] = -1;
        }
    }
}

/**
 * @brief    Aggregate membership degrees for fuzzy PID parameters
 * @details  Mathematical representation:
 *            μ_out[k] = Σ_iΣ_j (μ_e[i] × μ_ec[j])   for all rules where output = k
 *             -μ_e[i]: Membership degree of error in fuzzy set i
 *             -μ_ec[j]: Membership degree of error-change in fuzzy set j
 *             -μ_out[k]: Aggregated membership degree for output fuzzy set k
 * @return   void
 */
void FuzzyPID::GetSumGrad()
{
    // 清空隶属度累加数组(使用7个模糊集->NB到PB，故循环7次)
    for (int i = 0; i <= num_area - 1; i++)
    {
        KpgradSums[i] = 0;
        KigradSums[i] = 0;
        KdgradSums[i] = 0;

    }
    // 因为erro与erro_c有两个有效隶属度(左隶属度与右隶属度)，故总共有2×2种组合
    for (int i = 0; i < 2; i++)
    {
        if (e_grad_index[i] == -1)
        {
            continue;
        }
        for (int j = 0; j < 2; j++)
        {
            // 使用乘积推理法对隶属度进行合成(规则强度 = 左模糊集隶属度 * 右模糊集隶属度)
            if (ec_grad_index[j] != -1)
            {
                // 论域隶属度与实际索引的关系转换需要加3(NB=-3,而NB的索引为-3+3=0)
                int indexKp = Kp_rule_list[e_grad_index[i]][ec_grad_index[j]] + 3;
                int indexKi = Ki_rule_list[e_grad_index[i]][ec_grad_index[j]] + 3;
                int indexKd = Kd_rule_list[e_grad_index[i]][ec_grad_index[j]] + 3;

                // 自适应参数累加相同的数值，但是累加到不同的模糊集合
                //gradSums[index] = gradSums[index] + (e_gradmembership[i] * ec_gradmembership[j])* Kp_rule_list[e_grad_index[i]][ec_grad_index[j]];
                KpgradSums[indexKp] += e_gradmembership[i] * ec_gradmembership[j];
                KigradSums[indexKi] += e_gradmembership[i] * ec_gradmembership[j];
                KdgradSums[indexKd] += e_gradmembership[i] * ec_gradmembership[j];
            }
            else
            {
                continue;
            }

        }
    }

}

/**
 * @brief    Defuzzify aggregated membership degrees to crisp outputs
 * @details  The Center of Gravity formula for discrete sets is:
 *            y = [Σ(μ_i × v_i)] / Σ(μ_i)
 *             -μ_i: Membership degree of output fuzzy set i
 *             -v_i: Representative value (peak position) of fuzzy set i
 *             -y:   Crisp output value
 * @return void
 */
void FuzzyPID::GetOUT()
{
    // 利用重心法(Σ[隶属度×论域值])求解模糊论域值(由于反量化只影响缩放比例,不影响相对关系,故可省略分母)
    for (int i = 0; i < num_area - 1; i++)
    {
        qdetail_kp += kp_membership_values[i] * KpgradSums[i];
        qdetail_ki += ki_membership_values[i] * KigradSums[i];
        qdetail_kd+= kd_membership_values[i] * KdgradSums[i];
    }
}

/**
 * @brief      Execute one complete fuzzy PID control cycle
 * @details    This is the main entry point for the fuzzy PID controller.
 * @param[in]  e_max        Maximum expected error(physical units)
 * @param[in]  e_min        Minimum expected error(physical units)
 * @param[in]  ec_max       Maximum expected error change rate
 * @param[in]  ec_min       Minimum expected error change rate
 * @param[in]  kp_max       Maximum allowed Kp adjustment
 * @param[in]  kp_min       Minimum allowed Kp adjustment
 * @param[in]  erro         Current system error: e(t) = setpoint - measured
 * @param[in]  erro_c       Error change rate: de/dt ≈ (e(t) - e(t-1))/Δt
 * @param[in]  ki_max       Maximum allowed Ki adjustment
 * @param[in]  ki_min       Minimum allowed Ki adjustment
 * @param[in]  kd_max       Maximum allowed Kd adjustment
 * @param[in]  kd_min       Minimum allowed Kd adjustment
 * @param[in]  erro_pre     Error at previous time step: e(t-1)
 * @param[in]  erro_ppre    Error at time t-2: e(t-2)
 * @return     float        Control output(u)
 */
float FuzzyPID::FuzzyPIDcontroller(float e_max, float e_min, float ec_max, float ec_min, float kp_max, float kp_min, float erro, float erro_c,float ki_max,float ki_min,float kd_max,float kd_min,float erro_pre,float errp_ppre)
{
    errosum += erro;

    // 区间映射到模糊论域
    qerro = Quantization(e_max, e_min, erro);
    qerro_c = Quantization(ec_max, ec_min, erro_c);

    // 计算输入对各个论域的隶属度
    Get_grad_membership(qerro, qerro_c);

    // 模糊推理
    GetSumGrad();

    // 解模糊处理
    GetOUT();

    // 反量化(模糊论域 -> 实际值)
    detail_kp = Inverse_quantization(kp_max, kp_min, qdetail_kp);
    detail_ki = Inverse_quantization(ki_max, ki_min, qdetail_ki);
    detail_kd = Inverse_quantization(kd_max, kd_min, qdetail_kd);

    qdetail_kp = qdetail_ki = qdetail_kd = 0;

    if (qdetail_kp >= kp_max)
        qdetail_kp = kp_max;
    else if (qdetail_kp <= kp_min)
        qdetail_kp = kp_min;
    if (qdetail_ki >= ki_max)
        qdetail_ki = ki_max;
    else if (qdetail_ki <= ki_min)
        qdetail_ki = ki_min;
    if (qdetail_kd >= kd_max)
        qdetail_kd = kd_max;
    else if (qdetail_kd <= kd_min)
        qdetail_kd = kd_min;

    kp += detail_kp;
    ki += detail_ki;
    kd += detail_kd;

    // 根据极性,调整参数正负,具体问题具体分析
    if (kp < 0) kp = 0;
    if (ki < 0) ki = 0;
    if (kd < 0) kd = 0;

    detail_kp = detail_ki = detail_kd = 0;

    float output = kp * (erro - erro_pre) + ki * erro + kd * (erro - 2 * erro_pre + errp_ppre);
    return output;
}

/**
 * @brief      Map a value from physical range to normalized fuzzy domain
 * @details    The mapping is a linear transformation:
 *              y = [(x - minimum) / (max - min)] × 6 -3
 *               -x: Input value in physical units
 *               -y: Output value in normalized fuzzy domain [-3, 3]
 *               -minimum: Lower bound of physical range
 *               -maximum: Upper bound of physical range
 * @param[in]  maximum  Upper bound of physical input range
 * @param[in]  minimum  Lower bound of physical input range
 * @param[in]  x        Input value to quantize (in physical units)
 * @return     float  Quantized value in normalized range [-3, 3]
 */
float FuzzyPID::Quantization(float maximum,float minimum,float x)
{
    // 先归一化到[0, 1]区间,再扩展到[0, 6]区间,最后平移到[-3, 3]区间,实现区间映射
    float qvalues = 6.0f * (x - minimum) / (maximum - minimum) - 3;
    return qvalues;
}

/**
 * @brief      Inverse quantization: map from fuzzy domain to physical range
 * @details    The transformation is the inverse linear mapping:
 *              x = minimum + [(y + 3) / 6] × (maximum - minimum)
 *               -y: Input value in normalized fuzzy domain [-3, 3]
 *               -x: Output value in physical range [minimum, maximum]
 *               -minimum: Lower bound of physical range
 *               -maximum: Upper bound of physical range
 * @param[in]  maximum  Upper bound of target physical range
 * @param[in]  minimum  Lower bound of target physical range
 * @param[in]  y        Value in normalized fuzzy domain [-3, 3]
 * @return     float  Value mapped to physical range [minimum, maximum]
 */
float FuzzyPID::Inverse_quantization(float maximum, float minimum, float y)
{
    // [-3, 3]区间先平移到[0, 6]区间,再归一化到[0, 1]区间,扩展到实际物理范围,最后再平移到实际基准(即+minimum)实现反区间映射
    float x = (maximum - minimum) * (y + 3) / 6 + minimum;
    return x;
}