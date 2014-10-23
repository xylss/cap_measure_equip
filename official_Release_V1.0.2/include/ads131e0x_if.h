

/*	
 *	author	:		xuyl
 *	date	:		14-06-11
 *	filename:		ads131e0x_interface.h
 *	copyright:		xuyl 2014,6 --> 2019.6
 *  description:	it's designed for ADS131E08 spi interface operation
 */

 #ifndef __ADS131E0X_INTERFACE_H__
 #define __ADS131E0X_INTERFACE_H__
 

 /*ads131e0x's register*/
#define ADS131E0X_REG_ID_REG		0x00		//0xd2  ads131e08 ID
#define ADS131E0X_REG_CONFIG1		0x01
#define ADS131E0X_REG_CONFIG2		0x02
#define ADS131E0X_REG_CONFIG3		0x03
#define ADS131E0X_REG_FAULT			0x04

#define ADS131E0X_REG_CH1SET		0x05
#define ADS131E0X_REG_CH2SET		0x06
#define ADS131E0X_REG_CH3SET		0x07
#define ADS131E0X_REG_CH4SET		0x08
#define ADS131E0X_REG_CH5SET		0x09
#define ADS131E0X_REG_CH6SET		0x0A
#define ADS131E0X_REG_CH7SET		0x0B
#define ADS131E0X_REG_CH8SET		0x0C

#define ADS131E0X_REG_FAULT_STATP	0x12
#define ADS131E0X_REG_FAULT_STATN	0x13

#define  ADS131E0X_REG_GPIO			0x14


/*ads131e0x cmd define*/
//System Commands
#define ADS131SE0X_CMD_WAKEUP		0x02
#define ADS131E0X_CMD_STANDBY		0x04
#define ADS131E0X_CMD_RESET			0x06
#define ADS131E0X_CMD_START			0x08
#define ADS131E0X_CMD_STOP			0x0A
#define ADS131E0X_CMD_OFFSETCAL		0x1A
//Data Read Commands
#define ADS131E0X_CMD_RDATAC		0x10
#define ADS131E0X_CMD_SDATAC		0x11
#define ADS131E0X_CMD_RDATA			0x12
//Register Read Commands
#define ADS131E0X_CMD_RREG			0x20
#define ADS131E0X_CMD_WREG			0x40


#define ADS131E0X_TEST_NORMAL_INPUT 0x10
#define ADS131E0X_TEST_NOISE		0x11
#define ADS131E0X_TEST_MVDD			0x13
#define ADS131E0X_TEST_TEMPERATURE	0x14
#define ADS131E0X_TEST_SIG_DATA		0x15


// 小端模式, 需要先设置对齐,再反向, 默认32位对齐
typedef union{
	unsigned char spi_rd_buf[27];
	struct 
	{
		unsigned ads_state	: 24;
		unsigned channel_1	: 24;
		unsigned channel_2	: 24;
		unsigned channel_3	: 24;
		unsigned channel_4	: 24;
		unsigned channel_5	: 24;
		unsigned channel_6	: 24;
		unsigned channel_7	: 24;
		unsigned channel_8	: 24;	
	}ads_rlt;

}ad_val_t;
//

typedef struct{
	int ste;
	int ch_1;
	int ch_2;
	int ch_3;
	int ch_4;
	int ch_5;
	int ch_6;
	int ch_7;
	int ch_8;
}ad_cvs_rlt; //cvs: converse

int ads131e0x_send_SDATAC_cmd(void);

/*read 8 channel ad value to special array*/
int ads131e0x_read_ad_rlt(ad_cvs_rlt *prlt);

/*power on init, config some register, read AD value*/
int ads131e0x_power_init(void);
 

#endif

