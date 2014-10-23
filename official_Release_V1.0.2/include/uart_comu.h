


#ifndef __UART_COMU_H__
#define __UART_COMU_H__

#define NO_ACK				0
#define ACK_1				1
#define ACK_3				2


int uart_open(int ComPort, int baudrate, int databit,const char *stopbit, char parity);
void uart_close(int ComPort); 

int uart_write(int ComPort, char * data, int datalength);
int uart_read(int ComPort,void *data, int  datalength, int option);
int uart_clear(int ComPort);

#endif
