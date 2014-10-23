#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "table.h"
#include "demarcate.h"
#include "measure.h"

#define OVER_UPLIMIT_VALUE 2.2

typedef struct 
{
	int index;
	float ratio;
	float max_ratio;
}area_cmp_t;

float g_overflow = 0.0;

float adjust_follow_table(float cur_input, float *table, int len, float judge_range, float filter)
{
	float cur_output = cur_input;
	float tmp=0;
	int i;
	for(i=0; i<len; i++) {
		tmp = fabs(cur_input - table[i]);
		if(tmp < judge_range) {
			if(tmp > filter) {
				cur_output = table[i];
				break;
			}
		}			
	}
	
	return cur_output;
}

//单点进行区域筛选
int har_area_select(int times, float ori_ratio, parray_har_dem_t parray_dem) 
{
	int target_ratio = 0;
	int index = (times -3)/2;
	int i;
	
	for(i =0;i<40;i++) {
		if((ori_ratio > parray_dem[index][i].mid_low_range) && (ori_ratio < parray_dem[index][i].mid_high_range) ) {
			target_ratio = i+1;
			break;
		}	
	}
	
	return target_ratio;
}


int get_max_index(int times, area_cmp_t *parray, int num) 
{
	int i;
	int index = -1;
	int jump_flag = 0;
	int small_25_f = 0;
	area_cmp_t max_s;
	
	
	max_s.ratio = parray[0].ratio;
	max_s.index = parray[0].index;
	
	if(num == 0) 
		return -1;
	
	for(i=0;i<num;i++) {	
		/* jump_flag = 0; */
		/*filter lower ratio be detected situation, max_ratio must be small than index*/
		if(parray[i].index +1 < parray[i].max_ratio )
			continue;		
		
		if(max_s.ratio < parray[i].ratio) {		
			max_s.ratio = parray[i].ratio;
			max_s.index = parray[i].index;
			index = i;
		}	
	}
	
	if(max_s.ratio == parray[0].ratio) {
		index = 0;
	}
	
	return index ;
}


int har_section_area_compare(int times, float ori_mid, float ori_high, parray_har_dem_t parray_dem) 
{
	int target_ratio = 0;
	int index = (times -3)/2;
	int no = 0;
	int i,tmp;
	int cnt =0;
	area_cmp_t cmp_array[40];
	
	float overlap =0.0;
	float over_ratio =0.0;
	float total = 1.0;
	
	float reference_range = 0;
	float measure_range = ori_high - ori_mid;	
		
	memset(cmp_array, 0, sizeof(cmp_array));	

	for(i =0;i<40;i++) {
		reference_range = parray_dem[index][i].max_range - parray_dem[index][i].mid_high_range;
		total = (measure_range > reference_range)? measure_range:reference_range;
		
		//无区段交集情况
		if((ori_high < parray_dem[index][i].mid_high_range) || (ori_mid > parray_dem[index][i].max_range)) {
			continue;
		} 
		//小于参考区域 在其之间
		else if((ori_mid > parray_dem[index][i].mid_high_range) && (ori_high < parray_dem[index][i].max_range)) {
			overlap = ori_high - ori_mid;
			//total = parray_dem[index][i].max_range - parray_dem[index][i].mid_high_range;
			//over_ratio = 100 * overlap / total;
		} //参考区域被其包含
		else if((ori_mid < parray_dem[index][i].mid_high_range) && (ori_high > parray_dem[index][i].max_range)) {
			overlap = parray_dem[index][i].max_range - parray_dem[index][i].mid_high_range;
			//total = ori_high - ori_mid;
						
			over_ratio = 100 * overlap / total;
		}
		//在下限交叉 m左
		else if((ori_mid < parray_dem[index][i].mid_high_range) && (ori_high < parray_dem[index][i].max_range)) {
			overlap = ori_high - parray_dem[index][i].mid_high_range;
			total = ori_high - ori_mid;
			//total = parray_dem[index][i].max_range - parray_dem[index][i].mid_high_range;
			
		}//在上限交叉 m右
		else if((ori_mid > parray_dem[index][i].mid_high_range) && (ori_high > parray_dem[index][i].max_range)) {
			overlap = parray_dem[index][i].max_range - ori_mid;
			//total = parray_dem[index][i].max_range - parray_dem[index][i].mid_high_range;
			total = ori_high - ori_mid;
		}
		
		tmp = i+1;
		
		if(ori_high - g_overflow > tmp) continue;
		
		/*判断区段是否符合,大于90%做为判断依据*/
		over_ratio = 100 * overlap / total;		
				
		cmp_array[cnt].index = i;
		cmp_array[cnt].ratio = over_ratio;
		cmp_array[cnt].max_ratio = ori_high;
		cnt++;
		
		//printf(" cnt = %d, i =%d, over_ratio=%lf , overlap=%lf, total =%lf.\r\n", cnt, i, over_ratio, overlap, total);		
				
		if(over_ratio > 95) {
			/*过滤搜到到偏小区域(由于mid过小造成的), 要用max 进行筛选 max 应小于 target_ratio*/
			if(ori_high < tmp) {
				target_ratio = tmp;				
				break;
			} else {
				if(tmp %5 && ori_high - tmp < 2.5) {			
					target_ratio = tmp;			 	
					break;
				}
			}
		}
		
	}
	
	no = get_max_index(times, cmp_array, cnt);
	/*错误索引处理*/
	if(no == -1) {
		target_ratio = 0;
	}
	
	/*错误百分比处理*/
	if(cmp_array[no].ratio < 20) {
		target_ratio = 0;
	} else {
		target_ratio = cmp_array[no].index +1;
	}
	
	/*打印索引序号处理*/
	if(i == 40) {
		i = 41;
	}	
	
	//printf( i =%d, target_ratio = %d \r\n", i, target_ratio);	
	return target_ratio;
}

int harmonic_volt_calibration(phase_measure_t *phase, parray_har_dem_t parray)
{
	//printf("internal  谐波中文打印测试\r\n"); 
	int i,j;
	float *phar = NULL;
	float *pmax = NULL;
	float *pmid = NULL;	
	float mid_ratio = 0.0;
	float high_ratio = 0.0;
	float tmp_ratio =0;
	int index = 0;
	int rlt_ratio = 0.0;

	for(i=3; i<10; i+=2) {
		index = (i-3)/2;
		switch(i) {
		case 3:
			phar = &(phase->har_harmonic_3);
			pmax = &(phase->max_harmonic_3);
			pmid = &(phase->mid_harmonic_3);
			g_overflow = OVER_UPLIMIT_VALUE;
			break;
		case 5:
			phar = &(phase->har_harmonic_5);
			pmax = &(phase->max_harmonic_5);
			pmid = &(phase->mid_harmonic_5);
			g_overflow = 0;
			break;
		case 7:
			phar = &(phase->har_harmonic_7);
			pmax = &(phase->max_harmonic_7);
			pmid = &(phase->mid_harmonic_7);
			break;
		case 9:
			phar = &(phase->har_harmonic_9);
			pmax = &(phase->max_harmonic_9);
			pmid = &(phase->mid_harmonic_9);
			break;
		default:
			//printf("only deal with odd harmonic from 3 to 9, other ignore.\r\n");
			break;
		}
		
		/*judge is added harmonic by source*/
		if(*pmid < 0.15) 
			continue;	
		
		high_ratio = (100* *pmax / phase->mid_voltage_rms);
			
		mid_ratio = (100* *pmid / phase->mid_voltage_rms);
		
		//printf("%d times harmonic: *pmid =%lf V, mid_ratio = %lf, high_ratio =%lf. \r\n", i, *pmid, mid_ratio, high_ratio);
		
	#define HAR_CALI_PER_RATIO				
	#ifdef HAR_CALI_PER_RATIO		//以1%的谐波含量为刻度进行校准		
		// rlt_ratio = har_area_select(i, mid_ratio, parray);
		//区段比对, 90%重合即重合
		rlt_ratio = har_section_area_compare(i, mid_ratio, high_ratio, parray); 
		
		tmp_ratio = rlt_ratio;

		if(tmp_ratio >2.0) {
			rlt_ratio = (int)adjust_follow_table(tmp_ratio, har_tab, ARRAY_SIZE(har_tab), 2.1, 0.9);
		}

		//printf("after adjust_follow_table, rlt_ratio =%d. \r\n", rlt_ratio);
				
		if(rlt_ratio == 0) {			
			//printf("%d harmonics detected , voltage amplitude is between %lf\% ~ %lf\%, then calibrate to %d\% \r\n",i, parray[index][rlt_ratio].mid_high_range, parray[index][rlt_ratio].max_range, rlt_ratio);
		} else {	
			//printf("%d harmonics detected , voltage amplitude is between %lf\% ~ %lf\%, then calibrate to %d\% \r\n",i, parray[index][rlt_ratio -1].mid_high_range, parray[index][rlt_ratio -1].max_range, rlt_ratio);
		}
		 
		 *phar = 0.01* rlt_ratio * phase->mid_voltage_rms;
				
	#else 	// 5% 左右校准
		if(ratio > 1.5 && ratio <= 3) {
			printf("%d harmonics detected , voltage amplitude is between 1.5\% ~ 3\%, then calibrate to 3\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.03;
		} else if(ratio > 3 && ratio <= 5) {
			printf("%d harmonics detected , voltage amplitude is between 3\% ~ 5\%, then calibrate to 5\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.05;
		} else if(ratio > 5 && ratio <= 10) {
			printf("%d harmonics detected , voltage amplitude is between 5\% ~ 10\%, then calibrate to 10\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.1;
		} else if( ratio >10 && ratio <=15) {
			printf("%d harmonics detected , voltage amplitude is between 10\% ~ 15\%, then calibrate to 15\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.15;			
		} else if( ratio >15 && ratio <=20) {
			printf("%d harmonics detected , voltage amplitude is between 15\% ~ 20\%, then calibrate to 20\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.20;			
		} else if( ratio >20 && ratio <=25) {
			printf("%d harmonics detected , voltage amplitude is between 20\% ~ 25\%, then calibrate to 25\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.25;			
		} else if( ratio >25 && ratio <=30) {
			printf("%d harmonics detected , voltage amplitude is between 25\% ~ 30\%, then calibrate to 30\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.30;			
		} else if( ratio >30 && ratio <=35) {
			printf("%d harmonics detected , voltage amplitude is between 30\% ~ 35\%, then calibrate to 35\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.35;			
		} else if( ratio >35 && ratio <=40) {
			printf("%d harmonics detected , voltage amplitude is between 35\% ~ 40\%, then calibrate to 40\% \r\n",i);
			*pmax = phase->mid_voltage_rms * 0.40;			
		}
	#endif
		
	}
}