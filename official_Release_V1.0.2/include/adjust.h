
#ifndef __COMMON_H__
#define __COMMON_H__

void set_test(void);
float get_end_moment(void);
char get_phase_char(int i);
int quadrant_detect(float *pcur_buf, float *pvolt_buf, int len);
int converse_quadrant_angle(int quadrant, float *angle_buf, int len);

float numeric_sort(int num, float *pbuf );
float all_numeric_sort(int num, float *parray, float *pmin, float *pmax);
float calc_ave_sqrt(int* pData, int nNum);
float calculated_average(int len, float *pbuf);

float cur_adjust_follow_table(float cur_input, float *table, int len, float judge_range, float filter);
float freq_adjust_follow_table(float freq_input, float pre_freq, float *table, int len, float lvl_1, float lvl_2, float filter);
#endif