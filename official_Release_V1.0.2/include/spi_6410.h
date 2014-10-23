

#ifndef __SPI_TRANSFER_6410_H__
#define __SPI_TRANSFER_6410_H__

int appspi_send( char *buf, unsigned int len);
int appspi_recv_with_cmd( char cmd, char *tx_buf, int tx_len, char *rx_buf, int rx_len);
#endif

