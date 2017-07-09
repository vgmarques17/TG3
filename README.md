# TG3

Essa pasta contém arquivos relacionados ao projeto de TG3 sobre miografia por impedância elétrica, 
especialmente no que diz respeito ao controle da placa de aquisiçao


gcc -g -Wall -I. -o TG3_main TG3_main.c -L. Library/mylib.o -lmccusb  -L/usr/local/lib -lhidapi-libusb -lusb-1.0 -lm  -lgsl -lgslcblas -lconfig
