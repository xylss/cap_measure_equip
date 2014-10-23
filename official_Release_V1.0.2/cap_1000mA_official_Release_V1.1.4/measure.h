
#ifndef __MEASURE_H__
#define __MEASURE_H__

#define COUNT_MID_NUM  13	 		// 计算中值或平均值的数据
#define LVL2_COUNT_MID_NUM 9		// 二阶中值滤波计数

typedef struct 
{
	float K;
	float B;
}slope_t;

typedef struct 
{
	slope_t sec_9;
	slope_t sec_8;
	slope_t sec_7;
	slope_t sec_6;
	slope_t sec_5;
	slope_t sec_4;
	slope_t sec_3;
	slope_t sec_2;
	slope_t sec_1;
	slope_t sec_0;

}slope_section_t;	

typedef struct
{
	#define HAR_CALC_CYCLE 3  		 
	#define FREQ_FILTER_LEN 13
	#define FREQ_CALC_CYCLE 1		
	float volt_init_angle_buf[COUNT_MID_NUM+6];
	float cur_init_angle_buf[COUNT_MID_NUM+6];	
	
	float freq_buff[COUNT_MID_NUM * FREQ_CALC_CYCLE];
	float angle_buff[COUNT_MID_NUM];
	float vol_ave_sqrt_buf[COUNT_MID_NUM];
	float cur_ave_sqrt_buf[COUNT_MID_NUM];
	float volt_rms_buf[COUNT_MID_NUM];
	float cur_rms_buf[COUNT_MID_NUM];
	
	float harmonic_base_buf[COUNT_MID_NUM];
	float harmonic_3_buf[COUNT_MID_NUM];
	float harmonic_5_buf[COUNT_MID_NUM];
	float harmonic_7_buf[COUNT_MID_NUM];
	float harmonic_9_buf[COUNT_MID_NUM];
	float harmonic_base_angle_buf[COUNT_MID_NUM];

	float freq_filter_buff[FREQ_FILTER_LEN]; 			// 频率平滑处理缓冲数组
	float pre_freq; 									
	float ave_freq;
	float mid_freq;
	float lvl2_ave_freq;
	float lvl2_mid_freq;

	float pre_angle;
	float ave_angle;
	float mid_angle;
	float lvl2_ave_angle;
	float lvl2_mid_angle;

	float pre_volt_rms;
	float pre_cur_rms;
	float pre_volt_ave_sqrt;
	float pre_cur_ave_sqrt;
	
	float mid_harmonic_base;
	float mid_harmonic_3;	
	float mid_harmonic_5;
	float mid_harmonic_7;	
	float mid_harmonic_9;		
	float mid_harmonic_base_angle;
	
	float max_harmonic_3;
	float max_harmonic_5;
	float max_harmonic_7;
	float max_harmonic_9;	
	
	float pre_mid_harmonic_3;
	float pre_mid_harmonic_5;
	float pre_mid_harmonic_7;
	float pre_mid_harmonic_9;
	
	float pre_max_harmonic_3;
	float pre_max_harmonic_5;
	float pre_max_harmonic_7;
	float pre_max_harmonic_9;
	
	float har_harmonic_3;
	float har_harmonic_5;
	float har_harmonic_7;
	float har_harmonic_9;	
	float mid_volt_ave_sqrt;					// 均方根 对应的中值
	float mid_cur_ave_sqrt;
	float ave_volt_ave_sqrt;					// 均方根 对应的平均值
	float ave_cur_ave_sqrt;
	float lvl2_mid_volt_ave_sqrt;	
	float lvl2_mid_cur_ave_sqrt;	
	float lvl2_ave_volt_ave_sqrt;	
	float lvl2_ave_cur_ave_sqrt;	
	
	int quadrant;							 
	
	float dielectric_radian;					 
	float dielectric_angle;					 
	float dielectric_pf;						 
	float cos_angle;							 
	float cos_radian;							 
	float cos_pf;								 
	float capacity_current;					 
	float resistance_current;					 
	float voltage_rms;							 
	float current_rms;							 
	float capacitance;							
	
	float mid_voltage_rms;
	float mid_current_rms;
	float ave_voltage_rms;
	float ave_current_rms;
	
	unsigned int calculate_cnt;					// AD 计算结果计算
	unsigned int pre_calculate_cnt;
	unsigned int continue_cal_err_cnt;			//  调整	
	unsigned int lvl2_mid_calc_cnt; 			// 二阶中值滤波计数
	
	unsigned long angle_err_cnt;
	unsigned long freq_err_cnt;
	unsigned long volt_err_cnt;
	unsigned long cur_err_cnt ;
	unsigned long volt_rms_err_cnt;	
	unsigned long cur_rms_err_cnt;
	unsigned long total_filter_cnt;
		
	int filter_flag;		    				// 开启过滤标志 
	int err_cnt_flag;							// 用于错误计数赋值的标志, 置位后 才更新计数标准
	int phase_no_input_flag;
}phase_measure_t;
	
// 修正系数
typedef struct CORRECTION_COEFFICIENT
{
	float freq_coef;					
	float volt_rms_coef;				
	float cur_rms_coef;				
	float angle_coef;					
	float cos_angle_coef;				
	float pf_coef;						
	float cap_cur_coef;				
	float res_cur_coef;				
}corr_coeff_t;



#endif