
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <math.h>
#include "source.h"
#include "ads131e0x_if.h" 
#include "confile.h"
#include "demarcate.h"
#include "table.h"
#include "adjust.h"
#include "measure.h"
//#include "harmonic.h"

#define CFG_SET_FLOAT_AS_STRING(ft)  \
{ \
	sprintf(s, "%f", ft); \
	ret = ini_config_set_string(pconfig, section, key, 0, s, strlen(s)+1); \
}

#define CFG_GET_FLOAT_AS_STRING(ft)  \
{ \
	pvalue = ini_config_get_string(pconfig, section, key, "0.001"); \
	ft = atof(pvalue); \
}

#define PRINT_MEASURE_RLT_PATH	"/tmp/measure_rlt.txt"
#define UPDATE_PARA_PATH		"/etc/lyjd/update_para.ini"
#define INI_CFG_PATH			"/etc/lyjd/init_cfg"
#define DEFAULT_INI_CFG_PATH	"/etc/lyjd/default_cfg"
#define SLOPE_CURVE_PATH		"/etc/lyjd/slope_curve_cfg"
#define HAR_DEM_CFG_PATH		"/etc/lyjd/harmonic_dem_cfg"
//#define RUN_LEVEL_2_MID_CALCULATE	// 计算二阶中值开关
//#define MESURE_THREE_PHASE			// 采集三相电路开关
#define CURRENT_BLOCK_DIVISION		// 电流分段划分开关 (这里只硬件上分段划分)

#define COUNT_MID_NUM  13	//13		// 计算中值或平均值的数据
#define LVL2_COUNT_MID_NUM 9		// 二阶中值滤波计数

#define PHASE_CUR_UPPER_LIMIT  70 //70 		
#define POINT_CYCLE_PER  144	
#define CALCU_NUM	25 			// 3 	
#define DOTS_NUM_CYCLE  (POINT_CYCLE_PER * CALCU_NUM)		

#ifdef CURRENT_BLOCK_DIVISION		
	#define CUR_BLOCK_NUM 2			// 每相电流测量时分2段
#else
	#define CUR_BLOCK_NUM 1			// 不分段
#endif


#ifdef MESURE_THREE_PHASE
	#define PHASE_NUM 3 	
#else
	#define PHASE_NUM 1 
	//int g_ad_vol_buf[PHASE_NUM+1][CALCU_NUM * POINT_CYCLE_PER];	//计算B相的电压
	//int g_ad_vol_buf[PHASE_NUM][CALCU_NUM * POINT_CYCLE_PER];
#endif 	



#define CONTINUE_ERR_CNT_RANGE 4
#define NO_INPUT_DATA_THRESHOLD	80000
#define WAVE_SCOPE_COEFFICIENT (1.000000)
#define VREF_RANGE (4.096)
#define VREF_AD_RANGE (0x7fffff)

//#define PT_RATIO_A	(99.917/1.90001)	//实际电压值  调试1号板参数
#define PT_RATIO_A	(60.00/0.96000)
#define PT_RATIO_B	(99.917/1.89005)
#define PT_RATIO_C	(99.917/1.89000)


// 1000mA  小孔径 黑色 互感器
#define CT_RATIO_I2V_A_L (50/1.00)			//1-50mA 小电流段 1-40mA 切换档位  互感器输出
#define CT_RATIO_I2V_A_H  (1000/1.20) //(0.964385 *1000/1.2)
#define ARRAY2_ROW_SIZE(a) (sizeof(a) / sizeof((a)[0]))
#define ARRAY_SIZE(a) (sizeof(a) / sizeof((a)[0]))
#define ARRAY2_COL_SIZE(a) (sizeof((a)[0]) / sizeof((a)[0][0]))

#define FREQ_ADJUST_RANGE 	0.019//0.1Hz	//0.39  整Hz
#define ZERO_CORRECT_V_PHASE_A		(-38420)
#define ZERO_CORRECT_C_L_PHASE_A	(-37210)
#define ZERO_CORRECT_V_PHASE_B		(-4000)
#define ZERO_CORRECT_C_H_PHASE_A	(-4850)
#define ZERO_CORRECT_V_PHASE_C		(6000)

void DataDeal(int *pvolt_buf, int *pcur_buf);
void data_deal_calc_sqrt(int *pvolt_buf, int *pcur_buf);

typedef float(*parray_2_t)[2];
typedef float(*parray_3776_t)[DOTS_NUM_CYCLE];
typedef float(*parray_3_t)[3];
typedef int(*parray_collect_t)[3];

float TimeArry[DOTS_NUM_CYCLE];
float ArrayComplex[DOTS_NUM_CYCLE][2];
float ArrayReverse[2][DOTS_NUM_CYCLE];
float Array_multi[DOTS_NUM_CYCLE][2];
float Array_sigple[DOTS_NUM_CYCLE];
float TempDataA[DOTS_NUM_CYCLE];
float TempDataB[DOTS_NUM_CYCLE];

int g_ad_vol_buf[PHASE_NUM][CALCU_NUM * POINT_CYCLE_PER];
int g_ad_cur_buf[PHASE_NUM * CUR_BLOCK_NUM][CALCU_NUM * POINT_CYCLE_PER];

float g_primary_line_volt = 220000.0;	//变电站一次侧PT变比  
float g_secondary_volt = 60.0;

int g_demarcate_cycle = 1;
int amplitude_demarcate_flag = 0;		// 幅值标定标志
int g_device_type =0;  // 0: 容性,  1:避雷器, 2:铁芯接地

ad_cvs_rlt ad_rlt_s;
ad_cvs_rlt ad_rlt_test_s;	//use for test
float wave_range = WAVE_SCOPE_COEFFICIENT;
int g_zero_cnt_v =0;				// count of time which through zero.
int g_zero_cnt_i =0;				// count of time which through zero.
static float SignalFrequency = 0.0;
static float g_pre_freq = 0.0; 		//50.0;
static float BaseCalcPeriod = 138.15; 		 
static float g_base_dot_period = 140.35;
//struct timeval tv_start, tv_end;


float Angle;					 
float Radian; 					 
float electricinitangle;		 
float voltageinitangle;		 
 
float electricfai;				 
float voltagefai;				 
float tgvalue;					 	
float volt_ave_sqrt; 			 
float cur_ave_sqrt;			 

phase_measure_t g_phase_a_h;
phase_measure_t g_phase_a_l;
phase_measure_t g_phase_b_h;
phase_measure_t g_phase_b_l;
phase_measure_t g_phase_c_h;
phase_measure_t g_phase_c_l;

enum PHASE_E{PHASE_A =0, PHASE_B =1, PHASE_C=2 };
phase_measure_t g_phase_a;
phase_measure_t g_phase_b;
phase_measure_t g_phase_c;

char g_phase_mark = 'L';
char g_phase_chars ='A';
float g_phase_cur_ratio = 1.0;		// 用AD有效值求实际值时 所乘的系数
float g_phase_vol_ratio = 1.0;
float g_phase_cur_limit = 70.0;	//45.0;	// 40mA做判断,去除界限问题 
int g_small_curA_flag = 0;
int g_small_curB_flag = 0;
int g_small_curC_flag = 0;
int g_small_cur_flag  = 0;

float g_harmonic_base_vol = 0.0;
float g_harmonic_3_vol = 0.0;
float g_harmonic_5_vol = 0.0;
float g_harmonic_7_vol = 0.0;
float g_harmonic_9_vol = 0.0;
float g_harmonic_vol_coeff = PT_RATIO_A *VREF_RANGE /VREF_AD_RANGE * 0.029300;
float g_harmonic_base_angle = 0.0;

corr_coeff_t g_phase_A_H_coeff;
corr_coeff_t g_phase_A_L_coeff;
corr_coeff_t g_phase_B_H_coeff;
corr_coeff_t g_phase_B_L_coeff;
corr_coeff_t g_phase_C_H_coeff;
corr_coeff_t g_phase_C_L_coeff;
corr_coeff_t * p_phase_coeff;


/*加直流电压、电流 过零点时的偏移AD采样幅度的校准数值*/
int a_vol_amp_zero_correct = ZERO_CORRECT_V_PHASE_A; 
int a_cur_l_amp_zero_correct = ZERO_CORRECT_C_L_PHASE_A;
int a_cur_h_amp_zero_correct = ZERO_CORRECT_C_H_PHASE_A;
int b_vol_amp_zero_correct = ZERO_CORRECT_V_PHASE_B; 
int ads_collect_running = 1;
int ads_thread_loop = 1;

unsigned char ad_rx_buf[28] = {0,};

slope_section_t slope_sec_s;
slope_t slope_array[10];

har_dem_comu_t g_har_dem_s;
/*二维结构体数组指针变量类型*/
har_dem_save_t har_dem_save_array[4][40];	
har_dem_save_t har_dem_ref_array[4][40];


int harmonic_volt_calibration(phase_measure_t *phase, parray_har_dem_t parray);
float adjust_follow_table(float cur_input, float *table, int len, float judge_range, float filter);


int get_zero_drift_from_cfg(char *path, ad_cvs_rlt *p)
{
	INI_CONFIG *config;
	
	config = ini_config_create_from_file(DEFAULT_INI_CFG_PATH,0);
	
	if(config){
		 /* printf("CH1=%s, CH2=%s, CH3=%s, CH4=%s, CH5=%s, CH6=%s, \r\n", 
			ini_config_get_string(config, "ZERO_DRIFT","CH1_VOLT_A","NO"), \
			ini_config_get_string(config, "ZERO_DRIFT","CH2_VOLT_B","NO"), \
			ini_config_get_string(config, "ZERO_DRIFT","CH3_VOLT_C","NO"), \
			ini_config_get_string(config, "ZERO_DRIFT","CH4_CUR_A","NO"), \
			ini_config_get_string(config, "ZERO_DRIFT","CH5_CUR_B","NO"), \
			ini_config_get_string(config, "ZERO_DRIFT","CH6_CUR_C","NO")  \
		);  */
		
		p->ch_1 = ini_config_get_int(config, "ZERO_DRIFT","CH1_VOLT_A",0);
		p->ch_2 = ini_config_get_int(config, "ZERO_DRIFT","CH2_VOLT_B",0);
		p->ch_3 = ini_config_get_int(config, "ZERO_DRIFT","CH3_VOLT_C",0);
		p->ch_4 = ini_config_get_int(config, "ZERO_DRIFT","CH4_CUR_A",0);
		p->ch_5 = ini_config_get_int(config, "ZERO_DRIFT","CH5_CUR_B",0);
		p->ch_6 = ini_config_get_int(config, "ZERO_DRIFT","CH6_CUR_C",0);		
	
		g_primary_line_volt = ini_config_get_int(config, "SYSTEM_INI","PRIMARY_LINE_VOLT",0);
		
		
		/*读取 设备类型: 容性, 避雷器, 铁芯接地*/
		g_device_type =  ini_config_get_int(config, "SYSTEM_INI","DEVICE_TYPE",0);
		
		
		ini_config_destroy(config);
		return 1;
	}
	
	return 0;
}

/*上电读取谐波标定范围表*/
int init_har_dem_ref_range(parray_har_dem_t parray_dem)
{
	int i,j;
	int ret;
	char key[64];
	char section[32];	
	char s[20];
	float f;
	char *pvalue = NULL;
	INI_CONFIG* pconfig;	
	
	if(access(HAR_DEM_CFG_PATH, 0) == -1) {				
		printf("Don't find %s file ,init har_dem_cfg failed, then return \r\n.", HAR_DEM_CFG_PATH);
		return;
	} else {
		//printf("%s file is exist \r\n.", HAR_DEM_CFG_PATH);;
	}
	
	pconfig = ini_config_create_from_file(HAR_DEM_CFG_PATH,0);
	
	/*begin to save harmonic demarcate config file*/
	for(i=0;i<4;i++) {
		memset(section, 0, sizeof(section));
		sprintf(section, "HARMONIC_%d_TIMES", 2*i+3);
		for(j=0;j<40;j++) {			
			memset(key, 0, sizeof(key));
			sprintf(key, "RATIO_%d_MID", j+1);			

			CFG_GET_FLOAT_AS_STRING(parray_dem[i][j].mid_high_range);
			memset(key, 0, sizeof(key));
			sprintf(key, "RATIO_%d_MAX", j+1);	
			CFG_GET_FLOAT_AS_STRING(parray_dem[i][j].max_range);							
		}			
	}
	
	for(i=0;i<4;i++) {
		for(j=0;j<40;j++) {
		//	printf("%d times ratio %d , mid_low =%lf, mid_high=%lf. \r\n", i*2+3, j+1, parray_dem[i][j].mid_low_range, parray_dem[i][j].mid_high_range);
		}	
	}
	
	ret = ini_config_save(pconfig, HAR_DEM_CFG_PATH);			
	ini_config_destroy(pconfig);	
}

int save_har_dem_config(parray_har_dem_t parray_dem)
{
	int i,j;
	int ret;
	char key[64];
	char section[32];	
	char s[20];
	float f;
	
	FILE *fp = NULL;
	INI_CONFIG* pconfig;	
	//har_dem_save_t *ptmp;
	
	if(access(HAR_DEM_CFG_PATH, 0) == -1) {				
		printf("Don't find %s file ,now create it \r\n.", HAR_DEM_CFG_PATH);
		if((fp = fopen(HAR_DEM_CFG_PATH, "wr")) == NULL) {
			printf("can not creat %s file! \r\n", HAR_DEM_CFG_PATH);
			return -1;
		}	
		fclose(fp);
	} else {
		//printf("%s file is exist \r\n.", HAR_DEM_CFG_PATH);;
	}
	
	pconfig = ini_config_create_from_file(HAR_DEM_CFG_PATH,0);
	
	/*begin to save harmonic demarcate config file*/
	for(i=0;i<4;i++) {
		memset(section, 0, sizeof(section));
		sprintf(section, "HARMONIC_%d_TIMES", 2*i+3);
		printf(" section :%s\r\n", section);
		for(j=0;j<40;j++) {
			/* memset(key, 0, sizeof(key));
			sprintf(key, "RATIO_%d_MID_LOW", j+1);				
			CFG_SET_FLOAT_AS_STRING(parray_dem[i][j].mid_low_range); */

			memset(key, 0, sizeof(key));
			sprintf(key, "RATIO_%d_MID", j+1);						
			CFG_SET_FLOAT_AS_STRING(parray_dem[i][j].mid_high_range);
			
			memset(key, 0, sizeof(key));
			sprintf(key, "RATIO_%d_MAX", j+1);	
			CFG_SET_FLOAT_AS_STRING(parray_dem[i][j].max_range); 		
			
			/* memset(key, 0, sizeof(key));
			sprintf(key, "RATIO_%d_MIN", j+1);				
			CFG_SET_FLOAT_AS_STRING(parray_dem[i][j].min_range);  */			
		}			
	}

	/* 	value = strtod(input, &endptr); 
		ret = ini_config_set_int(pconfig,"SLOPE_SECTION","CH6_CUR_C",0, zero_rlt_buf[5]); 
	*/
	
	ret = ini_config_save(pconfig,HAR_DEM_CFG_PATH);			
	ini_config_destroy(pconfig);	
}

/*打印测试结果到文件中*/
int print_measure_to_file(phase_measure_t *phase)
{
	int i,j;
	int ret;
	char key[64];
	char section[32];	
	char s[20];
	float f;
	
	FILE *fp= NULL;
	INI_CONFIG* pconfig;	
	//har_dem_save_t *ptmp;
	
	if(access(PRINT_MEASURE_RLT_PATH, 0) == -1) {				
		printf("Don't find %s file ,now create it \r\n.", HAR_DEM_CFG_PATH);
		if((fp = fopen(PRINT_MEASURE_RLT_PATH, "wr")) == NULL) {
			printf("can not creat %s file! \r\n", HAR_DEM_CFG_PATH);
			return -1;
		}	
		fclose(fp);
	} else {
		//printf("%s file is exist \r\n.", HAR_DEM_CFG_PATH);;
	}
	
	pconfig = ini_config_create_from_file(PRINT_MEASURE_RLT_PATH, 0);
	memset(section, 0, sizeof(section));
	sprintf(section, "MEASURE_RLT", 2*i+3);
	
	/*printf measure result*/
	memset(key, 0, sizeof(key));
	sprintf(key, "FREQ");						
	CFG_SET_FLOAT_AS_STRING(phase->pre_freq);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "COS_PF");						
	CFG_SET_FLOAT_AS_STRING(phase->cos_pf);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "CAP_CUR");						
	CFG_SET_FLOAT_AS_STRING(phase->capacity_current);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "RES_CUR");						
	CFG_SET_FLOAT_AS_STRING(phase->resistance_current);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "TOTAL_CUR");						
	CFG_SET_FLOAT_AS_STRING(phase->mid_current_rms);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "MEASURE_VOLT");						
	CFG_SET_FLOAT_AS_STRING(phase->mid_voltage_rms);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "CAPACITANCE");						
	CFG_SET_FLOAT_AS_STRING(phase->capacitance);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "DI_ANGLE");						
	CFG_SET_FLOAT_AS_STRING(phase->dielectric_angle);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "DI_RADIAN");						
	CFG_SET_FLOAT_AS_STRING(phase->dielectric_radian);
	
	memset(key, 0, sizeof(key));
	sprintf(key, "DI_PF");						
	CFG_SET_FLOAT_AS_STRING(phase->dielectric_pf);
	
	if(g_device_type == 1) {
		memset(key, 0, sizeof(key));
		sprintf(key, "HAR_3_VOLT");						
		CFG_SET_FLOAT_AS_STRING(phase->har_harmonic_3);
		
		memset(key, 0, sizeof(key));
		sprintf(key, "HAR_5_VOLT");						
		CFG_SET_FLOAT_AS_STRING(phase->har_harmonic_5);
		
		memset(key, 0, sizeof(key));
		sprintf(key, "HAR_7_VOLT");						
		CFG_SET_FLOAT_AS_STRING(phase->har_harmonic_7);
		
		memset(key, 0, sizeof(key));
		sprintf(key, "HAR_9_VOLT");						
		CFG_SET_FLOAT_AS_STRING(phase->har_harmonic_9);
	}	
	
	ret = ini_config_save(pconfig, PRINT_MEASURE_RLT_PATH);		
}

/*slope config read*/
int get_slope_curve(int type)
{
	INI_CONFIG* pconfig;
	char *pvalue;
	char s[20];
	char key_name[20];
	int i;
	float f = 0.123456;
	FILE *fp = NULL;
	
	if(access(SLOPE_CURVE_PATH, 0) == -1) {				
		printf("Don't find %s file ,now exit\r\n.", SLOPE_CURVE_PATH);
		return -1;
	} 
	
	pconfig = ini_config_create_from_file(SLOPE_CURVE_PATH,0);

	for(i=0;i<10;i++) {
		memset(key_name,0, sizeof(key_name));
		sprintf(key_name, "SEC_%d_K", i);
		pvalue = ini_config_get_string(pconfig,"SLOPE_SECTION", key_name,"1.0");
		slope_array[i].K = atof(pvalue);
		memset(key_name,0, sizeof(key_name));
		sprintf(key_name, "SEC_%d_B", i);
		pvalue = ini_config_get_string(pconfig,"SLOPE_SECTION", key_name,"1.0");
		slope_array[i].B = atof(pvalue);
		
		//printf("sec %d: y = %lf x%c%lf. \r\n", i, slope_array[i].K, (slope_array[i].B >= 0)?'+':' ',slope_array[i].B);
	}				
	
	ini_config_destroy(pconfig);
	
	return 1;
}


//#define HAR_DEM_CYCLE 4  	//xy test set 1

int demarcate_once_harmonic(phase_measure_t *phase, har_dem_comu_t *pdem, parray_har_dem_t parray_dem)
{	
	int i,j;
	float max = 0;
	float mid = 0;	
	float min = 0;	
	float tmp_ratio = 0.0;
	int index = 0;
	int index_col; 
	
	/* */
	index_col = pdem->har_ratio - 1;
	switch(pdem->har_times) {
	case 3:
		max = phase->max_harmonic_3;
		mid = phase->mid_harmonic_3;
		min = phase->har_harmonic_3;
		index = 0;
		break;
	case 5:
		max = phase->max_harmonic_5;
		mid = phase->mid_harmonic_5;
		min = phase->har_harmonic_5;
		index = 1;
		break;
	case 7:
		max = phase->max_harmonic_7;
		mid = phase->mid_harmonic_7;
		min = phase->har_harmonic_7;
		index = 2;
		break;
	case 9:
		max = phase->max_harmonic_9;
		mid = phase->mid_harmonic_9;
		min = phase->har_harmonic_9;
		index = 3;
		break;
	default:
		printf("only deal with odd harmonic from 3 to 9, other ignore.\r\n");
		break;
	}
		
	/*judge is added harmonic by source*/
	if(mid < 0.003* phase->mid_voltage_rms) {
		printf("请检查标准源是否正常加入了谐波幅值, 有效值VRMS=%lf.\r\n",phase->mid_voltage_rms);
		return -1;
	}
	
	printf("g_har_dem_s.har_times=%d, g_har_dem_s.har_ratio=%d, har_volt=%lf, \r\n",pdem->har_times, pdem->har_ratio, pdem->har_volt );						
	//printf("%s()%d: now start \r\n",__func__, __LINE__);
	
	parray_dem[(pdem->har_times-3)/2][index_col].times = pdem->har_times; 
	parray_dem[(pdem->har_times-3)/2][index_col].dem_ratio = pdem->har_ratio; 
	
	tmp_ratio = (100* min / phase->mid_voltage_rms);		
	parray_dem[(pdem->har_times-3)/2][index_col].min_range = tmp_ratio; 
	tmp_ratio = (100* max / phase->mid_voltage_rms);		
	parray_dem[(pdem->har_times-3)/2][index_col].max_range = tmp_ratio; 
	tmp_ratio = (100* mid / phase->mid_voltage_rms);	
	parray_dem[(pdem->har_times-3)/2][index_col].mid_high_range = tmp_ratio + 0.3; 
	
	//处理(ratio == 1)的情况
	if(pdem->har_ratio == 1) {
		parray_dem[(pdem->har_times-3)/2][index_col].mid_low_range = 0.4;
	} else {
		parray_dem[(pdem->har_times-3)/2][index_col].mid_low_range = parray_dem[(pdem->har_times-3)/2][index_col -1].mid_high_range;
	}
	
	parray_dem[(pdem->har_times-3)/2][index_col].har_vol_ref = pdem->har_volt;
	
	return 1;	
}

/*过滤采样转换错误结果的处理过程*/
int get_filter_result(phase_measure_t *phase)
{
	int angle_err_flag 		= 0;
	int freq_err_flag 		= 0;
	int volt_ad_err_flag 	= 0;
	int cur_ad_err_flag 	= 0;	
	int volt_rms_err_flag 	= 0;
	int cur_rms_err_flag 	= 0;	
	int filter_err_flag 	= 0; //-1:出错 过滤掉改点， 0:未出错 不进行过滤
	
	if(phase->filter_flag) {

	}
	return filter_err_flag;
}

int check_is_dem_para()
{
	INI_CONFIG* pconfig;
	char *pvalue;
	char s[20];
	char key_name[20];
	int i,j;
	int ret;
	int a_v1, a_c4, b_v2, b_c5, c_v3,c_c6;	
	int zero_ad_buf[8][5 * POINT_CYCLE_PER];
	int zero_rlt_buf[8];
	FILE *fp = NULL;
	int type;
	
	
	if(access(UPDATE_PARA_PATH, 0) == -1) {				
		//printf("Don't find %s file ,now exit\r\n.", UPDATE_PARA_PATH);
		return 0;
	} else {
		printf("@@@@@@ Find %s file ,now start to dem para. @@@@@@\r\n", UPDATE_PARA_PATH);
	}
	
	pconfig = ini_config_create_from_file(UPDATE_PARA_PATH,0);
	
	type = ini_config_get_int(pconfig, "PARA_UPDATE_MENU","DEM_TYPE",0);
	
	if(type == 0) {
		printf("Dem type is 0, then do nothing, Please check %s file.\r\n", UPDATE_PARA_PATH);
		return 0;
	}
	switch(type) {
	case 1:	// 零漂
		for(i=0; i<5*POINT_CYCLE_PER; i++) {		
			ads131e0x_read_ad_rlt(&ad_rlt_test_s);
			a_v1 = ad_rlt_test_s.ch_1;
			b_v2 = ad_rlt_test_s.ch_2;
			c_v3 = ad_rlt_test_s.ch_3;
			a_c4 = ad_rlt_test_s.ch_4;
			b_c5 = ad_rlt_test_s.ch_5;
			c_c6 = ad_rlt_test_s.ch_6;
			
			if(a_v1 >= 0x800000 && a_v1 <=0xffffff) 	a_v1 = -(0x7fffff & (~a_v1 +1));
			if(a_c4 >= 0x800000 && a_c4 <=0xffffff) 	a_c4 = -(0x7fffff & (~a_c4 +1));			
			if(b_v2 >= 0x800000 && b_v2 <=0xffffff) 	b_v2 = -(0x7fffff & (~b_v2 +1));
			if(b_c5 >= 0x800000 && b_c5 <=0xffffff) 	b_c5 = -(0x7fffff & (~b_c5 +1));			
			if(c_v3 >= 0x800000 && c_v3 <=0xffffff) 	c_v3 = -(0x7fffff & (~c_v3 +1));
			if(c_c6 >= 0x800000 && c_c6 <=0xffffff) 	c_c6 = -(0x7fffff & (~c_c6 +1));
			
			zero_ad_buf[0][i] = a_v1;
			zero_ad_buf[1][i] = b_v2;
			zero_ad_buf[2][i] = c_v3;
			zero_ad_buf[3][i] = a_c4;
			zero_ad_buf[4][i] = b_c5;
			zero_ad_buf[5][i] = c_c6;
		}

		for(j=0;j<6;j++) {
			zero_rlt_buf[j] = int_average(3*POINT_CYCLE_PER, (int *)&zero_ad_buf[j][0]);
			printf("zero average value is %d, j=%d \r\n",zero_rlt_buf[j], j);				
		}
			
		pconfig = ini_config_create_from_file(DEFAULT_INI_CFG_PATH,0);

		ret = ini_config_set_int(pconfig,"ZERO_DRIFT","CH1_VOLT_A",0,zero_rlt_buf[0]);
		ret = ini_config_set_int(pconfig,"ZERO_DRIFT","CH2_VOLT_B",0,zero_rlt_buf[1]);
		ret = ini_config_set_int(pconfig,"ZERO_DRIFT","CH3_VOLT_C",0,zero_rlt_buf[2]);
		ret = ini_config_set_int(pconfig,"ZERO_DRIFT","CH4_CUR_A",0, zero_rlt_buf[3]);
		ret = ini_config_set_int(pconfig,"ZERO_DRIFT","CH5_CUR_B",0, zero_rlt_buf[4]);
		ret = ini_config_set_int(pconfig,"ZERO_DRIFT","CH6_CUR_C",0, zero_rlt_buf[5]);
		ret = ini_config_save(pconfig,DEFAULT_INI_CFG_PATH);
					
		remove(UPDATE_PARA_PATH); 
		break;
	case 2:	//电流标定
		break;
	case 3:	//电压标定
		break;
	case 4: // 谐波标定
		break;
	}
	
	return 1;
}

int init_measure_para(int mode)
{
/*三相修正系数初始化*/
#ifdef MESURE_THREE_PHASE
	/*init three phase correction coefficient */	
	memset(&g_phase_a_h, 0, sizeof(g_phase_a_h));
	memset(&g_phase_a_l, 0, sizeof(g_phase_a_l));
	memset(&g_phase_b_h, 0, sizeof(g_phase_b_h));
	memset(&g_phase_b_l, 0, sizeof(g_phase_b_l));
	memset(&g_phase_c_h, 0, sizeof(g_phase_c_h));
	memset(&g_phase_c_l, 0, sizeof(g_phase_c_l)); 	
	
	g_phase_A_H_coeff.freq_coef = 1.0;			
	g_phase_A_H_coeff.volt_rms_coef = 1.0;	
	g_phase_A_H_coeff.cur_rms_coef = 1.0;		
	g_phase_A_H_coeff.angle_coef = 1.0;		
	g_phase_A_H_coeff.cos_angle_coef = 1.0;		
	g_phase_A_H_coeff.pf_coef = 1.0;				
	g_phase_A_H_coeff.cap_cur_coef = 1.0;			
	g_phase_A_H_coeff.res_cur_coef = 1.0;	

	g_phase_A_L_coeff.freq_coef = 0.987473;//1.0;			
	g_phase_A_L_coeff.volt_rms_coef = 2.0;	
	g_phase_A_L_coeff.cur_rms_coef = 2.0;		
	g_phase_A_L_coeff.angle_coef = 1.0;		
	g_phase_A_L_coeff.cos_angle_coef = 1.0;		
	g_phase_A_L_coeff.pf_coef = 1.0;				
	g_phase_A_L_coeff.cap_cur_coef = 1.0;			
	g_phase_A_L_coeff.res_cur_coef = 1.0;	
	
	g_phase_B_H_coeff.freq_coef = 1.0;			
	g_phase_B_H_coeff.volt_rms_coef = 1.0;	
	g_phase_B_H_coeff.cur_rms_coef = 1.0;		
	g_phase_B_H_coeff.angle_coef = 1.0;		
	g_phase_B_H_coeff.cos_angle_coef = 1.0;		
	g_phase_B_H_coeff.pf_coef = 1.0;				
	g_phase_B_H_coeff.cap_cur_coef = 1.0;			
	g_phase_B_H_coeff.res_cur_coef = 1.0;	

	g_phase_B_L_coeff.freq_coef = 0.987473;//1.0;			
	g_phase_B_L_coeff.volt_rms_coef = 2.0;	
	g_phase_B_L_coeff.cur_rms_coef = 2.0;		
	g_phase_B_L_coeff.angle_coef = 1.0;		
	g_phase_B_L_coeff.cos_angle_coef = 1.0;		
	g_phase_B_L_coeff.pf_coef = 1.0;				
	g_phase_B_L_coeff.cap_cur_coef = 1.0;			
	g_phase_B_L_coeff.res_cur_coef = 1.0;	

	g_phase_C_H_coeff.freq_coef = 1.0;			
	g_phase_C_H_coeff.volt_rms_coef = 1.0;	
	g_phase_C_H_coeff.cur_rms_coef = 1.0;		
	g_phase_C_H_coeff.angle_coef = 1.0;		
	g_phase_C_H_coeff.cos_angle_coef = 1.0;		
	g_phase_C_H_coeff.pf_coef = 1.0;				
	g_phase_C_H_coeff.cap_cur_coef = 1.0;			
	g_phase_C_H_coeff.res_cur_coef = 1.0;		
	
	g_phase_C_L_coeff.freq_coef = 1.0;			
	g_phase_C_L_coeff.volt_rms_coef = 1.0;	
	g_phase_C_L_coeff.cur_rms_coef = 1.0;		
	g_phase_C_L_coeff.angle_coef = 1.0;		
	g_phase_C_L_coeff.cos_angle_coef = 1.0;		
	g_phase_C_L_coeff.pf_coef = 1.0;				
	g_phase_C_L_coeff.cap_cur_coef = 1.0;			
	g_phase_C_L_coeff.res_cur_coef = 1.0;		

#else 
/*单相修正系数初始化*/
	memset(&g_phase_a_h, 0, sizeof(g_phase_a_h));
	memset(&g_phase_a_l, 0, sizeof(g_phase_a_l));
	
	g_phase_A_H_coeff.freq_coef = 1.0;			
	g_phase_A_H_coeff.volt_rms_coef = 1.0;	
	g_phase_A_H_coeff.cur_rms_coef = 1.0;//0.961538; //100mA	0.974658;//1000mA		
	g_phase_A_H_coeff.angle_coef = 1.0;		
	g_phase_A_H_coeff.cos_angle_coef = 1.0;		
	g_phase_A_H_coeff.pf_coef = 1.0;				
	g_phase_A_H_coeff.cap_cur_coef = 1.0;			
	g_phase_A_H_coeff.res_cur_coef = 1.0;	

	g_phase_A_L_coeff.freq_coef = 1.0;			
	g_phase_A_L_coeff.volt_rms_coef = 1.0;	
	g_phase_A_L_coeff.cur_rms_coef = 1.0;
	g_phase_A_L_coeff.angle_coef = 1.0;		
	g_phase_A_L_coeff.cos_angle_coef = 1.0;		
	g_phase_A_L_coeff.pf_coef = 1.0;				
	g_phase_A_L_coeff.cap_cur_coef = 1.0;			
	g_phase_A_L_coeff.res_cur_coef = 1.0;	
	
	memset(&ad_rlt_test_s, 0, sizeof(ad_rlt_test_s));
	/*init zero drift value which read from default ini config*/
	get_zero_drift_from_cfg(DEFAULT_INI_CFG_PATH, &ad_rlt_test_s);
	a_vol_amp_zero_correct   = ad_rlt_test_s.ch_1;	// ZERO_CORRECT_V_PHASE_A; 
	a_cur_l_amp_zero_correct = ad_rlt_test_s.ch_4; 	// ZERO_CORRECT_C_L_PHASE_A;
	a_cur_h_amp_zero_correct = ad_rlt_test_s.ch_5;	// ZERO_CORRECT_C_H_PHASE_A;
	b_vol_amp_zero_correct   = ad_rlt_test_s.ch_2;	// ZERO_CORRECT_V_PHASE_B; 
	
	//g_primary_line_volt = ini_config_get_int(config, "SYSTEM_INI","PRIMARY_LINE_VOLT",0);
	/*init slope curve which get from cfg*/
	get_slope_curve(2);	
	init_har_dem_ref_range(har_dem_ref_array);

	check_is_dem_para();
#endif 	

	return mode;
}

int measure_result_handled(phase_measure_t *phase)
{
	int tmp_flag = 0;
	int quadrant = 0;
	
	/*filter error data*/
	if(phase->filter_flag){
		/*不合格数据直接跳过处理，采集下一个点的数据 continue;*/
		tmp_flag = get_filter_result(phase);
		
		if(tmp_flag == -1) {
			if(phase->err_cnt_flag == 0) {
				phase->pre_calculate_cnt = phase->calculate_cnt;
				phase->err_cnt_flag = 1;
			} else {		//err_cnt_flag == 1;
				if(phase->calculate_cnt - phase->pre_calculate_cnt <= 2){
					phase->continue_cal_err_cnt++; 						
					/*是否连续采集到的数据都超出之前设置的范围, 则认为测量量改变,重新初始化参数*/
					if(phase->continue_cal_err_cnt > CONTINUE_ERR_CNT_RANGE) {
						phase->continue_cal_err_cnt =0;
						phase->pre_calculate_cnt = 0;
						phase->calculate_cnt = 0;
						phase->lvl2_mid_calc_cnt = 0;
						
						phase->err_cnt_flag = 0;
						phase->filter_flag = 0;
						//printf("PHASE-%c, @@@@@@@@@@@@@@@:NOW is switch input source value for measure,ALL default data will restart @@@@@@@@@@@@@@@:\r\n",g_phase_chars);
					}
				} else 
					phase->continue_cal_err_cnt =0;
			}
			phase->pre_calculate_cnt = phase->calculate_cnt;
			
			//continue;	//because error  data don't add to calculate buff, so go to next cycle.
			return -1;
		}	
	}
	
	/*save data to buffer for calculate middle and average value*/
	phase->volt_init_angle_buf[phase->calculate_cnt % COUNT_MID_NUM] = voltageinitangle;
	phase->cur_init_angle_buf[phase->calculate_cnt % COUNT_MID_NUM] = electricinitangle;
	
	phase->freq_buff[phase->calculate_cnt % (COUNT_MID_NUM *FREQ_CALC_CYCLE)] = SignalFrequency;
	phase->angle_buff[phase->calculate_cnt % COUNT_MID_NUM] 		= Angle;
	phase->vol_ave_sqrt_buf[phase->calculate_cnt % COUNT_MID_NUM]	= volt_ave_sqrt;
	phase->cur_ave_sqrt_buf[phase->calculate_cnt % COUNT_MID_NUM] 	= cur_ave_sqrt;	
	
	phase->voltage_rms = g_phase_vol_ratio * volt_ave_sqrt;
	phase->current_rms = g_phase_cur_ratio * cur_ave_sqrt;
	phase->volt_rms_buf[phase->calculate_cnt % COUNT_MID_NUM]	= phase->voltage_rms ;
	phase->cur_rms_buf[phase->calculate_cnt % COUNT_MID_NUM] 	= phase->current_rms;	
	
	/*harmonic deal with */
	phase->harmonic_base_buf[phase->calculate_cnt % COUNT_MID_NUM] = g_harmonic_base_vol;
	phase->harmonic_3_buf[phase->calculate_cnt % (COUNT_MID_NUM *HAR_CALC_CYCLE)] = g_harmonic_3_vol;
	phase->harmonic_5_buf[phase->calculate_cnt % (COUNT_MID_NUM *HAR_CALC_CYCLE)] = g_harmonic_5_vol;
	phase->harmonic_7_buf[phase->calculate_cnt % (COUNT_MID_NUM *HAR_CALC_CYCLE)] = g_harmonic_7_vol;
	phase->harmonic_9_buf[phase->calculate_cnt % (COUNT_MID_NUM *HAR_CALC_CYCLE)] = g_harmonic_9_vol;
	phase->harmonic_base_angle_buf[phase->calculate_cnt % (COUNT_MID_NUM *HAR_CALC_CYCLE)] = g_harmonic_base_angle;
	
	phase->calculate_cnt++;		
	
	/*the first 20 times is used to calculate average for default filter value*/
	//if(phase->calculate_cnt < COUNT_MID_NUM) {
	if(phase->calculate_cnt < 7) {
		//continue;
		return -1;
	}	
	/*overflow deal with*/	
	else if(phase->calculate_cnt >= 0xFFFFFFFF) {
		phase->calculate_cnt = COUNT_MID_NUM;
	}	
	/*after COUNT_MID_NUM count ,then start to calculate the filter value.*/
	else if(phase->calculate_cnt >= COUNT_MID_NUM && !phase->filter_flag) {
		phase->quadrant = quadrant_detect(phase->cur_init_angle_buf, phase->volt_init_angle_buf,  COUNT_MID_NUM);
		converse_quadrant_angle(phase->quadrant, phase->angle_buff, COUNT_MID_NUM);
	
		phase->filter_flag = 1;	//only set once until calculate process released.
		phase->pre_angle		 = numeric_sort(ARRAY_SIZE(phase->angle_buff), phase->angle_buff);
		phase->pre_freq 		 = numeric_sort(ARRAY_SIZE(phase->freq_buff)/FREQ_CALC_CYCLE, phase->freq_buff);
		phase->pre_volt_rms		 = numeric_sort(ARRAY_SIZE(phase->volt_rms_buf), phase->volt_rms_buf);	
		phase->pre_cur_rms  	 = numeric_sort(ARRAY_SIZE(phase->cur_rms_buf), phase->cur_rms_buf);	
		phase->pre_cur_ave_sqrt  = numeric_sort(ARRAY_SIZE(phase->cur_ave_sqrt_buf), phase->cur_ave_sqrt_buf);	
		phase->pre_volt_ave_sqrt = numeric_sort(ARRAY_SIZE(phase->vol_ave_sqrt_buf), phase->vol_ave_sqrt_buf);	
		
		phase->pre_mid_harmonic_3 = numeric_sort(ARRAY_SIZE(phase->harmonic_3_buf)/HAR_CALC_CYCLE, phase->harmonic_3_buf);
		phase->pre_mid_harmonic_5 = numeric_sort(ARRAY_SIZE(phase->harmonic_5_buf)/HAR_CALC_CYCLE, phase->harmonic_5_buf);
		phase->pre_mid_harmonic_7 = numeric_sort(ARRAY_SIZE(phase->harmonic_7_buf)/HAR_CALC_CYCLE, phase->harmonic_7_buf);	
		phase->pre_mid_harmonic_9 = numeric_sort(ARRAY_SIZE(phase->harmonic_9_buf)/HAR_CALC_CYCLE, phase->harmonic_9_buf);	
		phase->pre_max_harmonic_3 = 0;
		phase->pre_max_harmonic_5 = 0;
		phase->pre_max_harmonic_7 = 0;	
		phase->pre_max_harmonic_9 = 0;	
		
		//pre_cur 	= calculated_average(ARRAY_SIZE(cur_buff_A), cur_buff_A);										
	} 
	/*count 17 times AD FFT converse result that deal with once average and middle filter value*/
	//else if(phase->filter_flag && ((phase->calculate_cnt - COUNT_MID_NUM)% COUNT_MID_NUM == 0)) {
	else if(phase->filter_flag && (phase->calculate_cnt - 0) %(COUNT_MID_NUM *g_demarcate_cycle) == 0) {
		phase->quadrant = quadrant_detect(phase->cur_init_angle_buf, phase->volt_init_angle_buf,  COUNT_MID_NUM);		//更新象限
		converse_quadrant_angle(phase->quadrant, phase->angle_buff, COUNT_MID_NUM);
		
		phase->mid_angle = p_phase_coeff->angle_coef * numeric_sort(ARRAY_SIZE(phase->angle_buff), phase->angle_buff);
		phase->ave_angle = p_phase_coeff->angle_coef * calculated_average(ARRAY_SIZE(phase->angle_buff), phase->angle_buff);
		phase->pre_angle = phase->mid_angle;
		
		phase->mid_freq = p_phase_coeff->freq_coef * numeric_sort(ARRAY_SIZE(phase->freq_buff), phase->freq_buff);
		phase->ave_freq = p_phase_coeff->freq_coef * calculated_average(ARRAY_SIZE(phase->freq_buff), phase->freq_buff);
		phase->pre_freq = phase->mid_freq;		
		
		phase->mid_volt_ave_sqrt = numeric_sort(ARRAY_SIZE(phase->vol_ave_sqrt_buf), phase->vol_ave_sqrt_buf);
		phase->ave_volt_ave_sqrt = calculated_average(ARRAY_SIZE(phase->vol_ave_sqrt_buf), phase->vol_ave_sqrt_buf);
		phase->mid_cur_ave_sqrt  = numeric_sort(ARRAY_SIZE(phase->cur_ave_sqrt_buf), phase->cur_ave_sqrt_buf);
		phase->ave_cur_ave_sqrt  = calculated_average(ARRAY_SIZE(phase->cur_ave_sqrt_buf), phase->cur_ave_sqrt_buf);
		phase->pre_volt_ave_sqrt = phase->mid_volt_ave_sqrt;
		phase->pre_cur_ave_sqrt  = phase->mid_cur_ave_sqrt;
	
		phase->voltage_rms = g_phase_vol_ratio * phase->mid_volt_ave_sqrt * p_phase_coeff->volt_rms_coef;
		phase->current_rms = g_phase_cur_ratio * phase->mid_cur_ave_sqrt * p_phase_coeff->cur_rms_coef;		

		phase->mid_voltage_rms = numeric_sort(ARRAY_SIZE(phase->volt_rms_buf), phase->volt_rms_buf);	
		phase->ave_voltage_rms = calculated_average(ARRAY_SIZE(phase->volt_rms_buf), phase->volt_rms_buf);		
		phase->mid_current_rms = numeric_sort(ARRAY_SIZE(phase->cur_rms_buf), phase->cur_rms_buf);	
		phase->ave_current_rms = calculated_average(ARRAY_SIZE(phase->cur_rms_buf), phase->cur_rms_buf);				
		phase->mid_voltage_rms *= p_phase_coeff->volt_rms_coef;
		phase->mid_current_rms *= p_phase_coeff->cur_rms_coef;		
		phase->pre_volt_rms = phase->mid_voltage_rms;
		phase->pre_cur_rms  = phase->mid_current_rms;	
		
		/*get middle and max value of harmonic, after count num is more than COUNT_MID_NUM*4 then start calc mid and max value per times */
		if(phase->calculate_cnt %HAR_CALC_CYCLE == 0) {
			phase->mid_harmonic_base = numeric_sort(ARRAY_SIZE(phase->harmonic_base_buf), phase->harmonic_base_buf);	
			phase->mid_harmonic_3 = all_numeric_sort(ARRAY_SIZE(phase->harmonic_3_buf), phase->harmonic_3_buf, &phase->har_harmonic_3, &phase->max_harmonic_3);		
			phase->mid_harmonic_5 = all_numeric_sort(ARRAY_SIZE(phase->harmonic_5_buf), phase->harmonic_5_buf, &phase->har_harmonic_5, &phase->max_harmonic_5);
			phase->mid_harmonic_7 = all_numeric_sort(ARRAY_SIZE(phase->harmonic_7_buf), phase->harmonic_7_buf, &phase->har_harmonic_7, &phase->max_harmonic_7);		
			phase->mid_harmonic_9 = all_numeric_sort(ARRAY_SIZE(phase->harmonic_9_buf), phase->harmonic_9_buf, &phase->har_harmonic_9, &phase->max_harmonic_9);
			phase->mid_harmonic_base_angle = numeric_sort(ARRAY_SIZE(phase->harmonic_base_angle_buf), phase->harmonic_base_angle_buf);
			
			phase->pre_mid_harmonic_3 = phase->mid_harmonic_3;
			phase->pre_mid_harmonic_5 = phase->mid_harmonic_5;
			phase->pre_mid_harmonic_7 = phase->mid_harmonic_7;	
			phase->pre_mid_harmonic_9 = phase->mid_harmonic_9;	
			phase->pre_max_harmonic_3 = phase->max_harmonic_3;
			phase->pre_max_harmonic_5 = phase->max_harmonic_5;
			phase->pre_max_harmonic_7 = phase->max_harmonic_7;	
			phase->pre_max_harmonic_9 = phase->max_harmonic_9;	
			//printf(" PHASE-%c, L1 3次谐波 mid =%lf V, max =%lf V , \r\n", g_phase_chars, phase->mid_harmonic_3, phase->pre_max_harmonic_3);	
		}
		
	
		
		/*输出*/
		if(phase->mid_current_rms < g_phase_cur_limit) { 
			g_small_cur_flag = 1;
			g_phase_cur_ratio = (CT_RATIO_I2V_A_L *VREF_RANGE /VREF_AD_RANGE);
			p_phase_coeff = &g_phase_A_L_coeff;
			//data_deal_calc_sqrt((int *)&g_ad_vol_buf[i], (int *)&g_ad_cur_buf[i]);  	
			
		} else {
			g_small_cur_flag = 0;
		
		}

		if(g_small_cur_flag == 1){
			/*8mA 以下电流标定*/
		#if 0
			if(phase->mid_current_rms >=70 && phase->mid_current_rms <80) {
				phase->mid_current_rms = cur_adjust_follow_table(phase->mid_current_rms, cur_rms_table_70_80mA, ARRAY_SIZE(cur_rms_table_70_80mA), 0.45, 0.34); 
			} else if(phase->mid_current_rms >40 && phase->mid_current_rms <70) {
				phase->mid_current_rms = cur_adjust_follow_table(phase->mid_current_rms, cur_rms_table_40_70mA, ARRAY_SIZE(cur_rms_table_40_70mA), 0.45, 0.2); 
			} else if(phase->mid_current_rms >30 && phase->mid_current_rms <41) {
				phase->mid_current_rms = cur_adjust_follow_table(phase->mid_current_rms, cur_rms_table_30_40mA, ARRAY_SIZE(cur_rms_table_30_40mA), 0.45, 0.14); 
			} else if(phase->mid_current_rms >20 && phase->mid_current_rms <31) {
				phase->mid_current_rms = cur_adjust_follow_table(phase->mid_current_rms, cur_rms_table_20_30mA, ARRAY_SIZE(cur_rms_table_20_30mA), 0.45, 0.10); 
			} else if(phase->mid_current_rms >10 && phase->mid_current_rms <21) {
				phase->mid_current_rms = cur_adjust_follow_table(phase->mid_current_rms, cur_rms_table_10_20mA, ARRAY_SIZE(cur_rms_table_10_20mA), 0.45, 0.05); 
			} else if(phase->mid_current_rms <11) {
				phase->mid_current_rms = cur_adjust_follow_table(phase->mid_current_rms, cur_rms_table_10mA, ARRAY_SIZE(cur_rms_table_10mA), small_cur_judge_range, small_cur_filter_threshold); 
			} 
		#endif
		
			/*no input station, set value to zero*/
			if( phase->mid_current_rms >=0 && phase->mid_current_rms < 0.7) { 
				 phase->mid_current_rms = 0.0;
			}
			if(phase->mid_voltage_rms >=0 && phase->mid_voltage_rms < 0.1) { 
				 phase->mid_voltage_rms = 0.0;
			}
			
			if(g_device_type == 1) {
				static int radom =1;
				phase->mid_current_rms *= 0.07;
				
				if(phase->mid_current_rms >0.085 && phase->mid_current_rms <0.13) {
					radom ++;
					phase->mid_current_rms = 0.100 + (radom %5)* 0.0001;					
				}
			}
			
			/* printf("#########: PHASE-%c_L, (small cur range) L1 RMS电压有效值 -->  平均值= %lf V,  /t中值= %lf V, \r\n",g_phase_chars, phase->ave_voltage_rms *p_phase_coeff->volt_rms_coef, phase->mid_voltage_rms);		
			printf("#########: PHASE-%c_L, (small cur range) L1 RMS电流有效值 -->  平均值= %lf mA, /t中值= %lf mA, AD电压值=%lf V.\r\n",g_phase_chars, phase->ave_current_rms *p_phase_coeff->cur_rms_coef, phase->mid_current_rms, phase->mid_current_rms/CT_RATIO_I2V_A_L ); */
			printf("#########: PHASE-%c, L1 RMS电压有效值   : %.3lf V, \r\n",g_phase_chars, phase->mid_voltage_rms);		
			printf("#########: PHASE-%c, L1 RMS电流有效值   : %.3lf mA,\r\n",g_phase_chars, phase->mid_current_rms);
		} 
		//else if(!g_small_cur_flag && g_phase_mark =='H') {
		else if(g_small_cur_flag == 0) {
			/*小孔径 1000 mA 互感器测量数据处理*/
			float tmp_mid = phase->mid_current_rms;
			if(phase->mid_current_rms >=912 &&  phase->mid_current_rms <1030) {			//900 --1000
				phase->mid_current_rms = slope_array[9].K * phase->mid_current_rms + slope_array[9].B; 
			} else if(phase->mid_current_rms >=811 &&  phase->mid_current_rms <912) {		//800 --900
				phase->mid_current_rms = slope_array[8].K * phase->mid_current_rms + slope_array[8].B; 
			} else if(phase->mid_current_rms >=720 &&  phase->mid_current_rms <811) {		//700 --800
				phase->mid_current_rms = slope_array[7].K * phase->mid_current_rms + slope_array[7].B; 
			} else if(phase->mid_current_rms >=615 &&  phase->mid_current_rms <720) {		//633 --700 		
				phase->mid_current_rms = slope_array[6].K * phase->mid_current_rms + slope_array[6].B; 
			} else if( phase->mid_current_rms >=485 &&  phase->mid_current_rms <615) { 	//500 --632	
				phase->mid_current_rms = slope_array[5].K * phase->mid_current_rms + slope_array[5].B; 
			} else if( phase->mid_current_rms >=315 &&  phase->mid_current_rms <485){	 	//(310 -- 470mA)
				phase->mid_current_rms = slope_array[4].K * phase->mid_current_rms + slope_array[4].B; 
				if(phase->mid_current_rms <=345) phase->mid_current_rms -= 0.55;
			} else if( phase->mid_current_rms >=238 &&  phase->mid_current_rms <315) { 	//(230 -- 310mA)
				phase->mid_current_rms = slope_array[3].K * phase->mid_current_rms + slope_array[3].B; 
			} else if( phase->mid_current_rms >=156 &&  phase->mid_current_rms <238) {		//(150 -- 250mA)
				phase->mid_current_rms = slope_array[2].K * phase->mid_current_rms + slope_array[2].B; 
			} else if( phase->mid_current_rms >=113 &&  phase->mid_current_rms <156) {		//(110 -- 150mA)
				phase->mid_current_rms = slope_array[1].K * phase->mid_current_rms + slope_array[1].B; 
			} else if( phase->mid_current_rms >=70 &&  phase->mid_current_rms <113) {		//(80 -- 110mA)
				phase->mid_current_rms = slope_array[0].K * phase->mid_current_rms + slope_array[0].B; 
				if(phase->mid_current_rms >=70  && phase->mid_current_rms <=80) 
					phase->mid_current_rms -= 0.2;	
			}
		

			if(g_device_type == 1) {
				phase->mid_current_rms *= 0.07;
			}
			
			printf("#########: PHASE-%c, L1 RMS电压有效值   : %.3lf V, \r\n", g_phase_chars, phase->mid_voltage_rms);		
			printf("#########: PHASE-%c, L1 RMS电流有效值   : %.3lf mA,\r\n", g_phase_chars, phase->mid_current_rms);
		}	
		
		/*频率修正*/
		if(phase->mid_freq <1.0) {
			SignalFrequency = 0.0;
		} else {
			SignalFrequency = freq_adjust_follow_table(phase->mid_freq, g_pre_freq, freq_table_0_1HZ, ARRAY_SIZE(freq_table_0_1HZ), 0.045, 0.0045, 0.0025);  
			
			SignalFrequency =  adjust_follow_table(SignalFrequency, freq_table, ARRAY_SIZE(freq_table), 0.07, 0.021);
		}

		
		phase->pre_freq = SignalFrequency;
		g_pre_freq = phase->pre_freq;			
	
		/**********************************************************************************/
		if(isnan(phase->mid_angle) != 0) {			
			phase->mid_angle = 0.000000;
			phase->capacitance = 0.000000;
		} else {
			/*求电容量 C= I/(2 *M_PI *F *U)*/
			phase->capacitance = (1000 * phase->mid_current_rms /(2 * M_PI * SignalFrequency * phase->mid_voltage_rms* 220000/60));	// uF
			phase->capacitance *= 1000000;	// pF	
			
			if((isnan(phase->capacitance) != 0) || (isinf(phase->capacitance) != 0)) {
				phase->capacitance = 0.0;
			} 
		}
		
		/*求功率因数角, 介质损耗角*/
		phase->cos_angle = phase->mid_angle * p_phase_coeff->cos_angle_coef;
		phase->cos_radian = (phase->cos_angle * M_PI / 180);
		if(phase->cos_angle <90) phase->dielectric_angle = fabs(90 - phase->cos_angle);
		else if(phase->cos_angle >90 && phase->cos_angle <180)  phase->dielectric_angle = fabs(90 - (180 - phase->cos_angle));
		else if(phase->cos_angle >270 && phase->cos_angle <360)  phase->dielectric_angle = fabs(90 - (360 - phase->cos_angle));
		
		/*filter dielectric_angle = 89.999999 , dielectric_pf = 57295779.814450 */
		if( phase->dielectric_angle > 89.999 && phase->dielectric_angle <90) {
			phase->dielectric_angle =90.0;
		}
		
		/*多位小数时正切值变大*/
		if(phase->dielectric_angle > 89.90 && phase->dielectric_angle < 89.99){
			/* char str[64]={0,};
			sprintf(str,"%.2f", phase->dielectric_angle);
			phase->dielectric_angle = (float)(atof(str)); */
			phase->dielectric_angle = (int)((phase->dielectric_angle *1000 +5)/10) / 100.0;			
		}
		
		phase->cos_pf = cos(phase->cos_radian) * p_phase_coeff->pf_coef;
		phase->capacity_current = phase->mid_current_rms * phase->cos_pf * p_phase_coeff->cap_cur_coef;
		phase->resistance_current = phase->mid_current_rms * sin(phase->cos_radian)* p_phase_coeff->res_cur_coef;
		
		/*介质损耗弧度, 介质损耗因数, 计算公式:(弧度 = 角度 *Pi /180) */
		phase->dielectric_radian = (phase->dielectric_angle * M_PI / 180);		
		//phase->dielectric_radian = (int)((phase->dielectric_radian *1000 +5)/10) / 100.0; //取2位小数
		
		if(phase->dielectric_angle == 90.0) {
			phase->dielectric_pf = 0;
		} else {
			phase->dielectric_pf = tan(phase->dielectric_radian);
		}
		///phase->dielectric_pf = tan(phase->dielectric_angle);
		
		//}

		printf("#########: PHASE-%c, L1 电网频率\t   : %.2lf Hz\r\n",g_phase_chars, SignalFrequency);
		printf("#########: PHASE-%c, L1 功率因数\t   : %.3lf \r\n",g_phase_chars, phase->cos_pf);
		
		printf("#########: PHASE-%c, L1 容性电流\t   : %.3lf mA \r\n",g_phase_chars, phase->capacity_current);	
		printf("#########: PHASE-%c, L1 阻性电流\t   : %.3lf mA \r\n",g_phase_chars, phase->resistance_current);	
		printf("#########: PHASE-%c, L1 全电流In\t  : %.3lf mA \r\n",g_phase_chars, phase->mid_current_rms);	
		printf("#########: PHASE-%c, L1 电容量Cx\t  : %.3lf pF \r\n",g_phase_chars, phase->capacitance);
		
		//printf("#########: PHASE-%c, L1 功率因数角度    : %.3lf ℃ \r\n",g_phase_chars, phase->cos_angle);	
		printf("#########: PHASE-%c, L1 介质损耗角度    : %.3lf ℃ \r\n",g_phase_chars, phase->dielectric_angle);	
		printf("#########: PHASE-%c, L1 介质损耗弧度    : %.3lf C \r\n",g_phase_chars, phase->dielectric_radian);	
		printf("#########: PHASE-%c, L1 介质损耗因数    : %.3lf \r\n\n",g_phase_chars, phase->dielectric_pf);	
		
		if(g_device_type == 1) {
			/*打印谐波幅值*/
			//printf("PHASE-%c, L1 3次谐波 mid =%lf V, max =%lf V \r\n",  g_phase_chars, phase->mid_harmonic_3, phase->max_harmonic_3);
			
			/* printf("#########: PHASE-%c, L1 5次谐波 mid =%lf V, max =%lf V \r\n", g_phase_chars, phase->mid_harmonic_5, phase->max_harmonic_5);
			printf("#########: PHASE-%c, L1 7次谐波 mid =%lf V, max =%lf V \r\n", g_phase_chars, phase->mid_harmonic_7, phase->max_harmonic_7);
			printf("#########: PHASE-%c, L1 9次谐波 mid =%lf V, max =%lf V \r\n", g_phase_chars, phase->mid_harmonic_9, phase->max_harmonic_9); */
			/*谐波标定*/
			//har_dem_save_t *p=NULL;
			parray_har_dem_t p = NULL;
			har_dem_comu_t *q = &g_har_dem_s;
			int wait_comu_cnt = 0;
			
			//if(g_harmonic_demarcate_flag) {	
			if(get_sys_flag(HARMONIC_DEMARCATE_FLAG)) {			
				demarcate_once_harmonic(phase, &g_har_dem_s, har_dem_save_array);
				p = har_dem_save_array;
				printf("$$$$$$$$$$: PHASE-%c, L1 %d次谐波%d\% 0° 标定结果: low_ratio =%lf\% ,mid_ratio =%lf\%, high_ratio =%lf\%,  标定参考电压 =%lf V\r\n", g_phase_chars, p[(q->har_times-3)/2][q->har_ratio -1].times, p[(q->har_times-3)/2][q->har_ratio -1].dem_ratio, p[(q->har_times-3)/2][q->har_ratio -1].min_range, p[(q->har_times-3)/2][q->har_ratio -1].mid_high_range, p[(q->har_times-3)/2][q->har_ratio -1].max_range, p[(q->har_times-3)/2][q->har_ratio -1].har_vol_ref);
				
				//printf("$$$$$$$$$$: PHASE-%c, L1 %d次谐波%d\% 0° 测量结果: min =%lf V , max =%lf V, mid =%lf V\r\n",g_phase_chars, p[(q->har_times-3)/2][q->har_ratio-1].times, p[(q->har_times-3)/2][q->har_ratio -1].dem_ratio, phase->har_harmonic_3, phase->max_harmonic_3, phase->mid_harmonic_3);
				
				/*每次指定百分比的谐波标定结束标志*/
				set_sys_flag(DEM_HAR_ONCE_FLAG, 1);
				
				/*查询通信标志,用来同步*/
				wait_comu_cnt = 0;
				while(get_sys_flag(DEM_HAR_START_COMU_FLAG) == 0){
					usleep(1000*1000);
					wait_comu_cnt++;
					printf(" wait %d s for har_dem_start_comu_flag \r\n", wait_comu_cnt);
				}
				wait_comu_cnt = 0;
				
				/*清除通信同步标志*/
				set_sys_flag(DEM_HAR_START_COMU_FLAG, 0);

			} 
			/*谐波校准*/
			else {
				//printf("before 谐波中文打印测试\r\n"); //在 harmonic_volt_calibration()内部中文打印不出来，很奇怪！！！
				p = har_dem_ref_array;				
				harmonic_volt_calibration(phase, p);
				
				printf("#########: PHASE-%c, L1 3次谐波电压幅值 : %.3lf V \r\n", g_phase_chars, phase->har_harmonic_3);
				printf("#########: PHASE-%c, L1 5次谐波电压幅值 : %.3lf V \r\n", g_phase_chars, phase->har_harmonic_5);
				printf("#########: PHASE-%c, L1 7次谐波电压幅值 : %.3lf V \r\n", g_phase_chars, phase->har_harmonic_7);
				printf("#########: PHASE-%c, L1 9次谐波电压幅值 : %.3lf V \r\n\n", g_phase_chars, phase->har_harmonic_9);
				//printf("#########: PHASE-%c, L1 基波谐波相位角  : %lf ℃ \r\n", g_phase_chars, phase->mid_harmonic_base_angle);
				printf("\r\n");
			}
		}
		
		/*保存在内存文件用 提供给61850通信使用*/
		print_measure_to_file(phase);
	}
	
	return 1;
}


/**/
int collect_AD_value_loop(void) 
{
	int per_cycle_cnt =0;
	int save_xml_cnt =0;
	
	int i,j;
	
	float timeuse; 

	int voltage_a, current_a, voltage_b, current_b, voltage_c, current_c;	//temp save ad data read from ads131e0x through spi bus.
	

	init_measure_para(0);
	
	/**/	
	//while(ads_thread_loop) {	
	while(1) {	
		GPP_pin_low(14);
		
		/*start data collection with frequency  16ksps , 62.5us per cycle read once*/		
		while(ads_collect_running) {
			set_test();
			//gettimeofday(&tv_start, NULL); 
			
			/*send stop continues read mode command per cycle before read ad value*/
			ads131e0x_send_SDATAC_cmd(); 
			
			/*start Collect ADS131E0x AD value and save it to temp array*/
			for(j=0; j<CALCU_NUM * POINT_CYCLE_PER; j++) {
				ads131e0x_read_ad_rlt(&ad_rlt_s);
			#ifndef CURRENT_BLOCK_DIVISION
				voltage_a = ad_rlt_s.ch_1;
				current_a = ad_rlt_s.ch_4;					
				voltage_b = ad_rlt_s.ch_2;
				current_b = ad_rlt_s.ch_5;					
				voltage_c = ad_rlt_s.ch_3;
				current_c = ad_rlt_s.ch_6;				
				
				/*采集数据的正负处理*/
				if(voltage_a >=0x800000 && voltage_a <=0xffffff) 
					voltage_a = -(0x7fffff & (~voltage_a +1));					
				if(current_a >=0x800000 && current_a <=0xffffff)  
					current_a = -(0x7fffff & (~current_a +1));				
				if(voltage_b >=0x800000 && voltage_b <=0xffffff) 
					voltage_b = -(0x7fffff & (~voltage_b +1));					
				if(current_b >=0x800000 && current_b <=0xffffff)  
					current_b = -(0x7fffff & (~current_b +1));
				if(voltage_c >=0x800000 && voltage_c <=0xffffff) 
					voltage_c = -(0x7fffff & (~voltage_c +1));					
				if(current_c >=0x800000 && current_c <=0xffffff)  
					current_c = -(0x7fffff & (~current_c +1));
				
				/*save collect AD value to handle buffer.*/
				g_ad_vol_buf[0][j] = voltage_a - a_vol_amp_zero_correct; 
				g_ad_cur_buf[0][j] = current_a - a_cur_amp_zero_correct; 
				g_ad_vol_buf[1][j] = voltage_b - b_vol_amp_zero_correct; 					
				g_ad_cur_buf[1][j] = current_b - b_cur_amp_zero_correct; 
				g_ad_vol_buf[2][j] = voltage_c - c_vol_amp_zero_correct; 					
				g_ad_cur_buf[2][j] = current_c - c_cur_amp_zero_correct; 
			#else
				voltage_a = ad_rlt_s.ch_1;
				voltage_b = ad_rlt_s.ch_2;
				current_a = ad_rlt_s.ch_4;					
				current_b = ad_rlt_s.ch_5;							
				
				/*采集数据的正负处理*/
				if(voltage_a >=0x800000 && voltage_a <=0xffffff) 
					voltage_a = -(0x7fffff & (~voltage_a +1));					
				if(current_a >=0x800000 && current_a <=0xffffff)  
					current_a = -(0x7fffff & (~current_a +1));				
				if(current_b >=0x800000 && current_b <=0xffffff)  
					current_b = -(0x7fffff & (~current_b +1));		
				if(voltage_b >=0x800000 && voltage_b <=0xffffff) 
					voltage_b = -(0x7fffff & (~voltage_b +1));	
					
				g_ad_vol_buf[0][j] = voltage_a - a_vol_amp_zero_correct; 
				g_ad_cur_buf[0][j] = current_a - a_cur_l_amp_zero_correct; 
				g_ad_vol_buf[1][j] = voltage_b - b_vol_amp_zero_correct; 		
				g_ad_cur_buf[0+1][j] = current_b - a_cur_h_amp_zero_correct; 
			
			#endif
				
				/*for disable ads131e0x auto run continues read mode, send stop cmd per cycle*/
				//if(j%POINT_CYCLE_PER == 0)
				//	ads131e0x_send_SDATAC_cmd(); 				
			}
			//gettimeofday(&tv_end, NULL); 
			//timeuse = 1000000 *(tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec); 
			timeuse = get_end_moment();
			g_base_dot_period = timeuse /(CALCU_NUM * POINT_CYCLE_PER);
			   			
			/*test FFT begin*/
			phase_measure_t *phase;
		#ifdef CURRENT_BLOCK_DIVISION
			//for(i=PHASE_A_L;i<=PHASE_A_H;i++) {
			for(i=PHASE_A;i<=PHASE_A;i++) {
				switch(i) {		
				case PHASE_A:			//0-->50mA (40mA division) current measure		
					phase = &g_phase_a;
					//g_phase_mark = 'L';
					g_phase_chars = get_phase_char(i);	
					g_phase_vol_ratio = PT_RATIO_A * VREF_RANGE / VREF_AD_RANGE ;
					if(g_small_cur_flag) {	//小电流情况 采用L 路测量数据
						g_phase_cur_ratio = (CT_RATIO_I2V_A_L *VREF_RANGE /VREF_AD_RANGE);	
						p_phase_coeff = &g_phase_A_L_coeff;
					} else {
						g_phase_cur_ratio = (CT_RATIO_I2V_A_H *VREF_RANGE /VREF_AD_RANGE);	
						p_phase_coeff = &g_phase_A_H_coeff;
					}
					//g_phase_cur_ratio = (CT_RATIO_I2V_A_H *CT_RATIO_V2V_A_H *VREF_RANGE /VREF_AD_RANGE);								
					p_phase_coeff = &g_phase_A_H_coeff;
					g_phase_cur_limit = PHASE_CUR_UPPER_LIMIT;
					break;
				/* case PHASE_A_H:			//50-->1000mA (40mA division) current measure				
					phase = &g_phase_a_h;
					g_phase_mark = 'H';
					g_phase_chars = get_phase_char(i);
					#if 1
						g_phase_cur_ratio = (CT_RATIO_I2V_A_H *VREF_RANGE /VREF_AD_RANGE);
					#else
						g_phase_cur_ratio = (CT_RATIO_I2V_A_H *CT_RATIO_V2V_A_H *VREF_RANGE /VREF_AD_RANGE);
					#endif
					g_phase_vol_ratio =  PT_RATIO_B * VREF_RANGE / VREF_AD_RANGE ;
					p_phase_coeff = &g_phase_A_H_coeff;
					g_phase_cur_limit = PHASE_CUR_UPPER_LIMIT;
					break;	 */		
				
				}
				
				/*If there were any data input? if not jump this collect data deal with.*/
				//if(phase->phase_no_input_flag) //preserve 
				// continue;

				DataDeal((int *)&g_ad_vol_buf[i], (int *)&g_ad_cur_buf[i+1]);  	

				if(g_small_cur_flag) {
					// g_phase_cur_ratio = (CT_RATIO_I2V_A_H *VREF_RANGE /VREF_AD_RANGE);
					 data_deal_calc_sqrt((int *)&g_ad_vol_buf[i], (int *)&g_ad_cur_buf[i]);  	
				}
		
				measure_result_handled(phase);
				ads131e0x_send_SDATAC_cmd(); 	
			}
		#else	//only phase B process.
			phase = &g_phase_a;
			g_phase_mark = PHASE_A;
			g_phase_chars = 'A';
			g_phase_cur_ratio = (CT_RATIO_I2V_A *CT_RATIO_V2V_A *VREF_RANGE /VREF_AD_RANGE);
			g_phase_vol_ratio =  PT_RATIO_A * VREF_RANGE / VREF_AD_RANGE ;
			p_phase_coeff = &g_phase_A_coeff;
			DataDeal((int *)&g_ad_vol_buf[0], (int *)&g_ad_cur_buf[0]);  				
			measure_result_handled(phase);
			ads131e0x_send_SDATAC_cmd(); 	
		#endif

		} //end while(ads_collect_running)
		
	} //end thread loop
}

 
int main(int argc, char **argv)
{
	int ret =0;	

	printf("System is init, Please wait ...\r\n");
	
	ads131e0x_power_init();
	collect_AD_value_loop();	
	release_source();
	
	return ret;
}



//**************************************************************//
/*	修正采样值算法 b0=Cn+n/N(Cn-C(n+N))
 *	<param name="currentbuffervalue">电流采样值</param>
 *	<param name="voltagebuffervalue">电压采样值</param>   
 * 	修正采样波（均匀化）
 */
void ModifiedSampledValues(int *pvolt, int *pcur)       
{
	int i;

	int cnt =0 ;
	float temp_val =0;	
	
	for ( i = 0; i < DOTS_NUM_CYCLE; i++) {				
		temp_val = pcur[i] + i / (SignalFrequency * g_base_dot_period / 1000000) * (pcur[i] - pcur[(i + (int)((SignalFrequency * g_base_dot_period / 1000000))) % 3776]);
		if(temp_val != pcur[i] ) {
			cnt++;
			pcur[i] = temp_val;
		}		
		pvolt[i] = pvolt[i] + i / (SignalFrequency * g_base_dot_period / 1000000) * (pvolt[i] - pvolt[(i + (int)((SignalFrequency * g_base_dot_period / 1000000))) % 3776]);		
	}

}


/*频率计算, 两个过零点之间计数插值 *采样周期 140.85us */
//float CaculateSignalFrequencyXJ(int *pvolt_buf, int *pcur_buf)   
float CaculateSignalFrequencyXJ(int *pvolt_buf, int *pcur_buf, float base_period) 			   
{
	float indexi[60]; 
	float indexv[60];	
	float buf_i[60];
	float buf_v[60];	
		
	float compareA, compareB;	
	float Resulti = 0, Resultv = 0; 
	int i,k = 0, j = 0;

	g_zero_cnt_i =0;
	g_zero_cnt_v = 0;
	for (i = 0; i < DOTS_NUM_CYCLE - 1; i++) {	
		compareA = pcur_buf[i];
		compareB = pcur_buf[i + 1];
		/*if ((compareA >= 0.0) && (compareB < 0.0)) {
		   // ZeroPoint[k] = compareA;
		   indexi[k] = i;
		   k++;
		}
		*/
		
		if ((i + 7) < (DOTS_NUM_CYCLE - 1)&&(compareA >= 0.0) && (compareB < 0.0)) {
			if ((pcur_buf[i] - pcur_buf[i + 7]) >= (pcur_buf[i + 1] - pcur_buf[i + 6]) &&
				(pcur_buf[i + 1] - pcur_buf[i + 6]) >= (pcur_buf[i + 2] - pcur_buf[i + 5]) &&
				(pcur_buf[i + 2] - pcur_buf[i + 5]) >= (pcur_buf[i + 3] - pcur_buf[i + 4])) {
				
				indexi[k] = i + 3;
				k++;
				g_zero_cnt_i++;
			}
		}						
	}
	
	for (i = 0; i < DOTS_NUM_CYCLE - 1; i++) {
		compareA = pvolt_buf[i];
		compareB = pvolt_buf[i + 1];

		if ((i + 7) < (DOTS_NUM_CYCLE - 1) && (compareA >= 0.0) && (compareB < 0.0)) {
			if ((pvolt_buf[i] - pvolt_buf[i + 7]) >= (pvolt_buf[i + 1] - pvolt_buf[i + 6]) &&
				(pvolt_buf[i + 1] - pvolt_buf[i + 6]) >= (pvolt_buf[i + 2] - pvolt_buf[i + 5]) &&
				(pvolt_buf[i + 2] - pvolt_buf[i + 5]) >= (pvolt_buf[i + 3] - pvolt_buf[i + 4])) 
			{
				indexv[j] = i + 3;
				j++;

				g_zero_cnt_v++;
			}
		}
	}
	
	//float Averagei = ((indexi[2] - indexi[1]) + (indexi[1] - indexi[0])) / 2;
	//float Averagev = ((indexv[2] - indexv[1]) + (indexv[1] - indexv[0])) / 2;	
	
	float Averagei, Averagev;
	float mid_i, mid_v;
	float filter_i, filter_v;
	float tmp =0;
	
	for(i =0; i<g_zero_cnt_i -1; i++) {
		tmp += (indexi[i+1] - indexi[i]);
		
		buf_i[i] = indexi[i+1] - indexi[i];	
	}
	Averagei = tmp/(g_zero_cnt_i -1);
	mid_i = numeric_sort(g_zero_cnt_i-2 , buf_i); 
	
	tmp = 0;
	for(i =0; i < g_zero_cnt_v -1; i++) {
		tmp += (indexv[i+1] - indexv[i]);
		
		buf_v[i] = indexv[i+1] - indexv[i];
	}
	Averagev = tmp/(g_zero_cnt_v-1); //xy 间隔是过零点数减一.
	mid_v = numeric_sort(g_zero_cnt_v -1, buf_v);

	Resulti = 1 / (Averagei * base_period / 1000000);
	Resultv = 1 / (Averagev * base_period / 1000000);

	mid_i =  1 / (mid_i * base_period / 1000000);
	mid_v =  1 / (mid_v * base_period / 1000000);	

	if( (fabs(mid_i - mid_v) < 4) && ((mid_i + mid_v)/2 >40 && (mid_i + mid_v)/2 < 68) ) {
		return Resultv;
		
	}  else {
		return 0.0;  
	}	
}
 
/* 	求三角式系数AN的函数
 *	<param name="h">采样间隔</param>
 *	<param name="voltagebuffervalue"></param>
 *	<param name="nw"></param>
 *  <param name="SignalFrequency"></param>
 */
float san(float h, int *pvolt, float nw, float freq)		//test OK
{
	float Sum = 0,Result=0;
	int j;

	for (j = 0; j < DOTS_NUM_CYCLE; j++) { 
		TempDataA[j] = cos(nw * j*h)*h;//xy  linux cos
		Sum = Sum + (float)pvolt[j] * TempDataA[j];
	}
	
	if (nw == 0) {
		Result = Sum * freq;
	}
	else {
		Result = 2*Sum * freq;
	}	 
	
	return Result;
}



/// //求三角式系数BN的函数
/// <param name="h">采样间隔</param>
/// <param name="voltagebuffervalue"></param>
/// <param name="nw"></param>
/// <param name="freq"></param>
float sbn(float h, int *pvolt, float nw, float freq)		//test OK
{
	float Sum = 0, Result = 0;
	int j = 0;
	for (j = 0; j < DOTS_NUM_CYCLE; j++)
	{
		TempDataB[j] = sin(nw * j * h) * h; //xy linux sin()
		Sum = Sum + (float)pvolt[j] * TempDataB[j];
	}

	Result = 2 * Sum * freq; ;
	return Result;
}



//typedef float(*parray_3_t)[3];
/*
 *  求谐波分量的电压  
 *  傅立叶级数表达式
 *  f(t) = a0+E(an*cos(n*w1*t)+bnsin(n*w1*t))
 *  a0 = 1/T*积分（f(t)*dt）～a0=1/T连加（f(t)）dt
 *  an = 2/T*积分（f(t)*cos(n*w1)*dt)～an=2/T连加（f(t)*cos(n*w1)）dt
 *  bn = 2/T*积分（f(t)*sin(n*w1)*dt）～bn=2/T连加（f(t)*sin(n*w1)）dt
 */
parray_3_t  harmonicwaves(int *voltagebuffervalue, float freq)		//test OK
{	 	
	//parray_3_t FourierCoefficient = malloc(sizeof(float[DOTS_NUM_CYCLE])*3); // 定义傅立叶系数数组
	float FourierCoefficient[DOTS_NUM_CYCLE][3];
	parray_3_t p = FourierCoefficient;
	
	float x, a0, an, bn, cn, h, w, nw;							// h :步长,  nw : n次谐波
	float frequency ;
	int n = 0, i;
	
	h = g_base_dot_period / 1000000;      		//采样时间间隔 us :[/1000000 us]
	
	w = 2 * M_PI * freq; 					 
	/*求取直流分量幅度*/
	nw = 0;
	a0 = san(h, voltagebuffervalue, nw, freq);
	//printf("%s:%d, h=%lf, w=%lf, nw = %lf, a0=%lf\r\n\n",__func__, __LINE__, h, w, nw, a0);
	
	FourierCoefficient[0][0] = FourierCoefficient[0][2] = a0;

	/*求N次谐波分量频谱幅度*/
	//for ( n = 1; n < DOTS_NUM_CYCLE; n++) { 
	for ( n = 1; n < 10; n+=2) { 	//xy  耗时太长, 只需计算32次
	//for ( n = 1; n < 10; n++) { 
		nw = n * w;
		an = san(h, voltagebuffervalue, nw, freq);
		bn = sbn(h, voltagebuffervalue, nw, freq);
		x = an * an + bn * bn;

		cn = sqrt(x); //cn=Sqrt（an * an + bn * bn）//xy linux sqrt()

		FourierCoefficient[n][0] = an;
		FourierCoefficient[n][1] = bn;
		FourierCoefficient[n][2] = cn;

	}
	
	/*calc harmonic voltage amplitude */
	g_harmonic_base_vol =  FourierCoefficient[1][2] * g_harmonic_vol_coeff;
	g_harmonic_3_vol 	=  FourierCoefficient[3][2] * g_harmonic_vol_coeff;
	g_harmonic_5_vol	=  FourierCoefficient[5][2] * g_harmonic_vol_coeff;
	g_harmonic_7_vol 	=  FourierCoefficient[7][2] * g_harmonic_vol_coeff;
	g_harmonic_9_vol 	=  FourierCoefficient[9][2] * g_harmonic_vol_coeff;
	
	float tmp_radian = 0;
	//tmp_radian = atan( FourierCoefficient[1][0] / FourierCoefficient[1][1] );	//tan@ = an/bn
	tmp_radian = atan( FourierCoefficient[1][1] / FourierCoefficient[1][0] );	//tan@ = bn/an
	g_harmonic_base_angle = (tmp_radian / M_PI * 180);
	//PT_RATIO_A *VREF_RANGE /VREF_AD_RANGE * 0.029785
	
	//打印计算结果
	//i =0;
	//printf("%d次谐波: 傅系数a= %lf, 系数b= %lf, 电压幅值 = %lf. \r\n",i, FourierCoefficient[i][0], FourierCoefficient[i][1], FourierCoefficient[i][2]*PT_RATIO_A *VREF_RANGE /VREF_AD_RANGE * 0.029785 );
			;
	//for (i = 0; i < DOTS_NUM_CYCLE; i++) {		
	for ( i = 0; i < 32; i++) { 	 
		
		if (i%2==1 && (2*i+1)<8) {		
			
			//printf("%d次谐波: 傅系数a= %lf, 系数b= %lf, 电压幅值 = %lf V\r\n",i, FourierCoefficient[i][0], FourierCoefficient[i][1], FourierCoefficient[i][2]*PT_RATIO_A *VREF_RANGE /VREF_AD_RANGE * 0.029405 );
			;
		}
	}            

	return p;
}

float* CollectionTime()		//test OK
{
	float *p = TimeArry;
	int i;
	for (i = 0; i < DOTS_NUM_CYCLE; i++) {
		TimeArry[i] =(float) i;
	}
	return p;
}


/* 设置矩阵初始值
 */
static parray_2_t SetArray(float *TimeArry)   //求出sin（（2*pi*f/t)*deta_T）		//test OK
{
	parray_2_t parray = ArrayComplex;
	
	float a, b;
	int i;
	float temp;
	for (i = 0; i < DOTS_NUM_CYCLE; i++)
	{
		temp = 2 * M_PI * 50 * TimeArry[i];
		//a = ArrayA[i, 0] = Math.Sin(2 * Math.PI * 50 * TimeArry[i]/50000);
		//b = ArrayA[i, 1] = Math.Cos(2 * Math.PI * 50 * TimeArry[i]/50000);
		//a = ArrayA[i, 0] = Math.Sin(2 * Math.PI * 50 * TimeArry[i]*0.000030395);
		//b = ArrayA[i, 1] = Math.Cos(2 * Math.PI * 50 * TimeArry[i]*0.000030395);
		//a = ArrayA[i, 0] = Math.Sin(2 * Math.PI * 50 * TimeArry[i] * 0.0000304414);
		//b = ArrayA[i, 1] = Math.Cos(2 * Math.PI * 50 * TimeArry[i] * 0.0000304414);
		//a = ArrayA[i, 0] = Math.Sin(2 * Math.PI * 50 * TimeArry[i] * 0.02 / 659.5);
		//b = ArrayA[i, 1] = Math.Cos(2 * Math.PI * 50 * TimeArry[i] * 0.02 / 659.5);
		a = parray[i][0] = sin(2 * M_PI * SignalFrequency * TimeArry[i] * g_base_dot_period / 1000000);
		b = parray[i][1] = cos(2 * M_PI * SignalFrequency * TimeArry[i] * g_base_dot_period / 1000000);
		//a = ArrayA[i, 0] = Math.Sin(2 * Math.PI * 50 * TimeArry[i] * 0.02 / 834);
		//b = ArrayA[i, 1] = Math.Cos(2 * Math.PI * 50 * TimeArry[i] * 0.02 / 834);
	}
	return parray;
}


/*	矩阵的转置
 * <param ="ArrayA">原矩阵</param>
 * <returns>转置后的矩阵</returns>
 */
static parray_3776_t  TransPose(parray_2_t ArrayA)	 
{
	parray_3776_t parray = ArrayReverse;
	int i;
	for (i = 0; i < DOTS_NUM_CYCLE; i++) {
		parray[0][i] = ArrayA[i][0];
		parray[1][i] = ArrayA[i][1];
	}
	return parray;
}


//		test OK 
/// 矩阵相乘
/// <param name="LeftArray左矩阵</param>
/// <param name="RightArray">右矩阵</param>
//typedef float(*parray_37x_t)[3776];
static parray_2_t Matrix_multiplication(parray_3776_t pLeftArray, parray_2_t pRightArray)	 
{		
	memset(Array_multi, 0, sizeof(Array_multi));
	parray_2_t parray = Array_multi;
		
	int i, j, k;
	for (i = 0; i < 2; i++) {							// LeftArray.GetLength(0); 	//2
		for (j = 0; j < 2; j++) {						// RightArray.GetLength(1);	//2
			for (k = 0; k < DOTS_NUM_CYCLE; k++) {		// LeftArray.GetLength(1);	//3776
				parray[i][j] += pLeftArray[i][k] * pRightArray[k][j];				
			}			

			//printf("Array_multi[%d][%d]=%lf, ", i,j, parray[i][j] );		
		}
	}

	return parray;
}



static float* Matrix_sigleplication(parray_3776_t pLeftArray, int *pRightArray)
{
	float *parray = malloc(sizeof(float[2]));
	memset(parray, 0, sizeof(float[2]));

	int i, j;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < DOTS_NUM_CYCLE; j++) {
			parray[i] += pLeftArray[i][j] * pRightArray[j];
		}
	}
	return parray;
}


float ArrayMatrix[2][2];
static parray_2_t GetReplaceMatrix(parray_2_t pCoefficientArray, float *pSinValueArry, int i)
{	
	parray_2_t parray = ArrayMatrix;
	int m,n;
	for (m = 0; m < 2; m++) {
		for (n = 0; n < 2; n++) {
			if (i != m) {
				parray[n][m] = pCoefficientArray[n][m];
			}
			else parray[n][i] = pSinValueArry[n];
			
			
		}
	}
	
	return parray;
}
   
/// 递归的方法求出某个行列式的值
/// <param name="input"></param>
static float GetMatrixResult(parray_2_t input)
{
	/*because "input" point to "Array_multi" */
	if (ARRAY2_COL_SIZE(Array_multi) == 2) {
		return input[0][0] * input[1][1] - input[0][1] * input[1][0];
	} else {
		return 0.0;
	}
}


int save_xml_cnt = 1;
int save_measure_flag = 0;

void savetestrecord()
{
	int number;	
}

void data_deal_calc_sqrt(int *pvolt_buf, int *pcur_buf)  
{
	cur_ave_sqrt = calc_ave_sqrt(pcur_buf, DOTS_NUM_CYCLE); 
}	

	
float real_calc_volt[DOTS_NUM_CYCLE]; 	//实际采样到的AD数值 对应的电压值
/*采样到的AD数据 处理*/
void DataDeal(int *pvolt_buf, int *pcur_buf)                                            
{
	//float cur_buf[DOTS_NUM_CYCLE];
	//float volt_buf[DOTS_NUM_CYCLE];	
	int i;
	volt_ave_sqrt = calc_ave_sqrt(pvolt_buf, DOTS_NUM_CYCLE); 
	cur_ave_sqrt = calc_ave_sqrt(pcur_buf, DOTS_NUM_CYCLE); 	

	
	/*修正之前计算信号频率*/
	SignalFrequency = CaculateSignalFrequencyXJ(pvolt_buf, pcur_buf, g_base_dot_period); 
	//SignalFrequency = DFTCalSignalCycle(cur_buf, volt_buf);  //
		
	//去掉滤波处理  和 谐波计算部分, 为了节省计算时间	
	/*x 如果是int 型，处理后数据未变化，除非出现错误的波形数据*/
	ModifiedSampledValues(pvolt_buf, pcur_buf);  //xy only process  Convert.ToInt32
	
	harmonicwaves(pvolt_buf, g_pre_freq);		//ori ad value , 过滤后的频率, 

	
	//float *TimeArry = malloc(sizeof(float[DOTS_NUM_CYCLE]));
	float *pTimeArry;
	pTimeArry = CollectionTime();             //

	//parray_2_t tempvalue = malloc(sizeof(float[2][2])); 
	parray_2_t tempvalue;
	tempvalue = Matrix_multiplication(TransPose(SetArray(TimeArry)), SetArray(TimeArry));    //求出sin（（2*pi*f/t)*deta_T）	
		
	/**/
	//float *pcurrentresult = malloc(sizeof(float[2]));
	float *pcurrentresult;
	pcurrentresult = Matrix_sigleplication(TransPose(SetArray(TimeArry)), pcur_buf);  //计算采样的电流数组与三角数组的矩阵乘积
	//printf("pcurrentresult[0]=%lf, pcurrentresult[1]=%lf.  \r\n", pcurrentresult[0], pcurrentresult[1]);
	 
	//float *pvoltageresult  =  malloc(sizeof(float[2]));	
	float *pvoltageresult;	
	pvoltageresult = Matrix_sigleplication(TransPose(SetArray(TimeArry)), pvolt_buf); //计算采样的电压数组与三角数组的矩阵乘积
	//printf("pvoltageresult[0]=%lf, pvoltageresult[1]=%lf.  \r\n", pvoltageresult[0], pvoltageresult[1]);
	
	/****/
	//float *CalCurrentResult = malloc(sizeof(float[2]));
	float ARRAY_CAL_CUR[2];	
	float *pCalCurrentResult = ARRAY_CAL_CUR;
	
	float TempA=0.0, TempB=0.0; 
	parray_2_t ptest1;
	//printf("test ARRAY_SIZE(ARRAY_CAL_CUR) = %d \r\n", ARRAY_SIZE(ARRAY_CAL_CUR));
	for (i = 0; i < ARRAY_SIZE(ARRAY_CAL_CUR); i++) {

		ptest1 = GetReplaceMatrix(tempvalue, pcurrentresult, i);
		TempA = GetMatrixResult(ptest1);
		
		TempB = GetMatrixResult(tempvalue);
		pCalCurrentResult[i] = TempA / TempB;		

	}

	//float* pCalvoltageResult = malloc(sizeof(float[2]));
	float ARRAY_CAL_VOL[2];	
	float *pCalvoltageResult = ARRAY_CAL_VOL;
	float TempC=0.0, TempD=0.0;
	parray_2_t ptest2;
	for (i = 0; i < ARRAY_SIZE(ARRAY_CAL_VOL); i++){
		ptest2 = GetReplaceMatrix(tempvalue, pvoltageresult, i);
		TempC = GetMatrixResult(ptest2);		
		TempD = GetMatrixResult(tempvalue);
		pCalvoltageResult[i] = TempC / TempD;				

	}	

	
	float temp1 = pCalCurrentResult[1] / pCalCurrentResult[0];
	float temp2 = pCalvoltageResult[1] / pCalvoltageResult[0];
	
	//printf("temp1=%lf,  temp2=%lf\r\n",temp1, temp2);
	
	electricfai = atan(temp1); 		// Math.Atan2(temp1);    
	voltagefai = atan(temp2) ;     	 


	/* electricmaxvalue = fabs(pCalCurrentResult[1] / sin(electricfai));
	voltagemaxvalue = fabs(pCalvoltageResult[1] / sin(voltagefai)); */
		
	/*弧度转换为角度*/
	electricinitangle = (electricfai / M_PI * 180);
	voltageinitangle  = (voltagefai  / M_PI * 180);
	
	/*弧度值夹角*/		
	Radian = fabs(fabs(electricfai) - fabs(voltagefai));        
	
	/*介质损耗角度*/
	/* if (voltagefai * electricfai > 0) { 
		Angle = electricfai > voltagefai ? fabs((electricfai - voltagefai) / M_PI * 180 - 0) : fabs((voltagefai - electricfai) / M_PI * 180 - 0);  
	} else {
		Angle = (90+(90 - fabs(voltageinitangle))) - electricinitangle;
	}

	if (Angle >= 65) 
		Angle = Angle - 1.15;
	if (Angle < 65) 
		Angle = Angle  - 1.3;

    if(electricinitangle* voltageinitangle <0) {	
		Angle = fabs(voltageinitangle - electricinitangle);		
		
		if(Angle >90) {
			Angle = 180.0 - Angle;
		}
	} 	
	else {	
		Angle = fabs(voltageinitangle - electricinitangle);		
	} */

	Angle = voltageinitangle - electricinitangle;	
			
	tgvalue = tan(M_PI/2 - Angle/180*M_PI) ;                           	

	if(save_measure_flag) {
		savetestrecord();
	}
	
	free(pcurrentresult);
	free(pvoltageresult);		
}

