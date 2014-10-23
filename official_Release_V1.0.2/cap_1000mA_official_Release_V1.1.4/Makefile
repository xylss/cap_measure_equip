
TARGET = cap_equip.exe
INSTALL_DIR = /exports

.PHONY:all clean $(TARGET)

$(TARGET): main.c 
	arm-linux-gcc -o $(TARGET)  main.c table.c harmonic.c -I../include -L../lib -lcfg -lgpio -lADS131E0x_official  -lpthread -lm -ladjust -luart -ldemarcate

	
#-lpthread 
	
install:
	cp $(TARGET) $(INSTALL_DIR)
	
clean:
	rm $(TARGET)
	
#-I/mnt/hgfs/vm_win8_share/spi-6410-transfer/
	

