


#ifndef __DEMARCATE_H__
#define __DEMARCATE_H__

#define STD_COMU_TYPE_DEMARCATE 1
#define STD_COMU_TYPE_MEASURE 	2

typedef struct
{
	//har :  harmonic
	int   har_ratio; 		// (1-->40)
	char  har_times;		// (2-->32)
	float har_volt;			
	float har_angle;		

}har_dem_comu_t;

typedef struct
{
	float max_range;
	float min_range;
	float mid_high_range;
	float mid_low_range;	
	
	int dem_ratio;		//ref, output result (1-->40)
	int times;
	float har_vol_ref;
}har_dem_save_t;

typedef har_dem_save_t(*parray_har_dem_t)[40];	

#define DEM_HAR_ONCE_FLAG 		2
#define DEM_HAR_START_COMU_FLAG	3
#define DEM_HAR_FINISH_FLAG 	4
#define HARMONIC_DEMARCATE_FLAG 6
int get_sys_flag(int type);
int set_sys_flag(int type , int value);

/* int is_litte_endian(void);
unsigned char xor_sum_check(char *buf, char len);
int cmd_parse(void);

int set_gear(char va, char vb, char vc, char ia, char ib, char ic);
int set_disp_gui(char type);
int set_amplitude_value(float va, float vb, float vc, float ia, float ib, float ic);
int amplitude_demarcate_comu(int type);

int set_harmonic_value(char channel, char count, char harmonic, float amplitude, float angle); */
int harmonic_demarcate_std_comu(har_dem_comu_t *pdem, int type);
int harmonic_demarcate_process(har_dem_comu_t *pdem, int type);
int clear_terminal(int index);


#endif
