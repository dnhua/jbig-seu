#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "jbig.h"
#include "jbig_table.h"


#define JBGMAXSIZE (5*1024*1024*10)
#define true 1
#define false 0





typedef int bool;
typedef char uint8_t;
typedef short int uint16_t;


size_t len1 = 4960;
unsigned char bitmap1[4500000] = { 0 };

int read2buffer(unsigned char *buffer, size_t size, size_t nmemb, unsigned char *jbigbuf, size_t *index)
{
	memcpy(buffer, jbigbuf + (*index)*size, size*nmemb);
	*index = *index + nmemb;
	return nmemb;
}

bool Jbigmodule(uint8_t* inbuf, uint8_t* outbuf, uint16_t width)
{

	int i, j, result;
	struct jbg_dec_state s;
	unsigned char *buffer, *p;
	size_t buflen, len, cnt;
	size_t datalen = 450000;
	size_t index = 0;
	unsigned long xmax = 4294967295UL, ymax = 4294967295UL, max;
	int plane = 1, use_graycode = 1, diagnose = 0, multi = 0;
	buflen = 8000;
	buffer = (unsigned char *)malloc(buflen);
	if (!buffer) {
		printf("Sorry, not enough memory available!\n");
		exit(1);
	}

	//intial
	jbg_dec_init(&s);
	jbg_dec_maxsize(&s, xmax, ymax);
	/* read BIH first to check VLENGTH */	
	len = read2buffer(buffer, 1, 20, inbuf, &index);
	if (len < 20) {
		exit(1);
	}
	result = JBG_EAGAIN;
	do {
		cnt = 0;
		p = (unsigned char *)buffer;	
		while (len > 0 &&
			(result == JBG_EAGAIN || (result == JBG_EOK && multi))) {
			result = jbg_dec_in(&s, p, len, &cnt);
			if (result != JBG_EOK && result != JBG_EAGAIN) {
				printf("JBIGDECODE_FAILED %d\n", result);
				return false;
			}
			p += cnt;
			len -= cnt;
		}
		if (!(result == JBG_EAGAIN || (result == JBG_EOK && multi))) {
			printf("JBIG_OK\n");
			break;
		}
		len = read2buffer(buffer, 1, buflen, inbuf, &index);
		datalen -= len;	
	} while (datalen > 0);

	//write image to pmb file	
	FILE *fout = NULL;
	fout = fopen("./image/0001.pbm", "wb");
	if (plane >= 0 && jbg_dec_getplanes(&s) <= plane) {
		printf("begin to write image\n");
		fprintf(fout, "P4\n%ld %ld\n", jbg_dec_getwidth(&s),
			jbg_dec_getheight(&s));
		if (fout == NULL) {
			printf("file is NULL\n");
		}
		else {
			printf("write image data\n");
			outbuf = jbg_dec_getimage(&s, 0);


			fwrite(outbuf, 1, jbg_dec_getsize(&s), fout);
			fclose(fout);
			printf("write image data done!\n");
		}
		jbg_dec_free(&s);
		free(buffer);
	}
	return true;
}

bool decodeJbig(char * infilename)
{


	unsigned char jbigbuf[450000] = { 0 };
	unsigned char *outbuf;
	outbuf = (unsigned char *)malloc(JBGMAXSIZE);
	unsigned char *addr;
	FILE *fp;
	int i = 0;
	if ((fp = fopen(infilename, "r")) == NULL)
	{
		printf("!!can not open the file.");
		return false;
	}
	addr = jbigbuf;
	printf("begin read image\n");
	while (!feof(fp))
	{
		*addr = fgetc(fp);
		addr++;
		i++;
	}
	fclose(fp);
	printf("read image end\n");

	printf("begin decode\n");
	bool flag = Jbigmodule(jbigbuf, outbuf, 8);


	printf("decode end \n");
	return flag;
}

int main()
{
	printf("---------- begin decode first image ----------\n");
	char infilename[50] = "./image/yzfjbg.jbg";
	bool flag = decodeJbig(infilename);
	printf("\n---------- first image end ----------\n");

	
	system("pause");
	return 0;
}

