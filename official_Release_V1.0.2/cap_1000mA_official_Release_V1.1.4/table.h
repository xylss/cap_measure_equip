

#ifndef __ATABLE_H__
#define __ATABLE_H__

//计算二维数组是 计算序号0 即多少ROW e.g, array[m][n] --> (ARRAY2_ROW_SIZE(a) == m)
#define ARRAY2_ROW_SIZE(a) (sizeof(a) / sizeof((a)[0]))
#define ARRAY_SIZE(a) (sizeof(a) / sizeof((a)[0]))
//计算二维数组是 计算序号1 即多少COL e.g, array[m][n] --> (ARRAY2_ROW_SIZE(a) == n)
#define ARRAY2_COL_SIZE(a) (sizeof((a)[0]) / sizeof((a)[0][0]))

extern float freq_table_0_1HZ[230];
extern float freq_table[24];
extern float har_tab[9];
#endif
