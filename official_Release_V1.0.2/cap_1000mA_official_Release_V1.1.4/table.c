



float small_cur_judge_range = 0.3; //0.3mA 
float small_cur_filter_threshold = 0.02;
float cur_rms_table_10mA[] = {
1,2,3,4,5,6,7,8,9,10
};
float cur_rms_table_10_20mA[] = {
11,12,13,14,15,16,17,18,19
};
float cur_rms_table_20_30mA[] = {
21,22,23,24,25,26,27,28,29
};
float cur_rms_table_30_40mA[] = {
31,32,33,34,35,36,37,38,39
};
float cur_rms_table_40_70mA[] = {
40,41,42,43,44,45,46,47,48,49,
50,51,52,53,54,55,56,57,58,59,
60,61,62,63,64,65,66,67,68,69
};
float cur_rms_table_70_80mA[] = {
70,71,72,73,74,75,76,77,78,79,
};

float freq_table[24] = {
43.00, 44.00, 45.00, 46.00, 47.00, 
48.00, 49.00, 50.00, 51.00, 52.00, 
53.00, 54.00, 55.00, 56.00, 57.00, 
58.00, 59.00, 60.00, 61.00, 62.00, 
63.00, 64.00, 65.00, 66.00, 
};
float har_tab[9] = {
0, 5, 10, 15, 20, 25, 30, 35, 40
};

float freq_table_0_1HZ[230] = {
44.00, 44.10, 44.20, 44.30, 44.40, 44.50, 44.60, 44.70, 44.80, 44.90,
45.00, 45.10, 45.20, 45.30, 45.40, 45.50, 45.60, 45.70, 45.80, 45.90, 
46.00, 46.10, 46.20, 46.30, 46.40, 46.50, 46.60, 46.70, 46.80, 46.90,
47.00, 47.10, 47.20, 47.30, 47.40, 47.50, 47.60, 47.70, 47.80, 47.90,
48.00, 48.10, 48.20, 48.30, 48.40, 48.50, 48.60, 48.70, 48.80, 48.90,
49.00, 49.10, 49.20, 49.30, 49.40, 49.50, 49.60, 49.70, 49.80, 49.90,
50.00, 50.10, 50.20, 50.30, 50.40, 50.50, 50.60, 50.70, 50.80, 50.90,
51.00, 51.10, 51.20, 51.30, 51.40, 51.50, 51.60, 51.70, 51.80, 51.90,
52.00, 52.10, 52.20, 52.30, 52.40, 52.50, 52.60, 52.70, 52.80, 52.90,
53.00, 53.10, 53.20, 53.30, 53.40, 53.50, 53.60, 53.70, 53.80, 53.90,
54.00, 54.10, 54.20, 54.30, 54.40, 54.50, 54.60, 54.70, 54.80, 54.90,
55.00, 55.10, 55.20, 55.30, 55.40, 55.50, 55.60, 55.70, 55.80, 55.90,
56.00, 56.10, 56.20, 56.30, 56.40, 56.50, 56.60, 56.70, 56.80, 56.90,
57.00, 57.10, 57.20, 57.30, 57.40, 57.50, 57.60, 57.70, 57.80, 57.90,

58.00, 58.10, 58.20, 58.30, 58.40, 58.50, 58.60, 58.70, 58.80, 58.90,
59.00, 59.10, 59.20, 59.30, 59.40, 59.50, 59.60, 59.70, 59.80, 59.90,
60.00, 60.10, 60.20, 60.30, 60.40, 60.50, 60.60, 60.70, 60.80, 60.90,
61.00, 61.10, 61.20, 61.30, 61.40, 61.50, 61.60, 61.70, 61.80, 61.90,
62.00, 62.10, 62.20, 62.30, 62.40, 62.50, 62.60, 62.70, 62.80, 62.90,
63.00, 63.10, 63.20, 63.30, 63.40, 63.50, 63.60, 63.70, 63.80, 63.90,
64.00, 64.10, 64.20, 64.30, 64.40, 64.50, 64.60, 64.70, 64.80, 64.90,
65.00, 65.10, 65.20, 65.30, 65.40, 65.50, 65.60, 65.70, 65.80, 65.90,
66.00, 66.10, 66.20, 66.30, 66.40, 66.50, 66.60, 66.70, 66.80, 66.90,
};

