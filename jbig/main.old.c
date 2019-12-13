//#include <stdio.h>
//#include "jbig.h"
//#include "jbig_table.h"
//
//void output_bie(unsigned char *start, size_t len, void *file)
//{
//	fwrite(start, 1, len, (FILE *)file);
//
//	return;
//}
//
//int main()
//{
//	printf("---------- start ----------\n");
//	unsigned char bitmap[15] = {
//		/* 23 x 5 pixels, "JBIG" */
//		0x7c, 0xe2, 0x38, 0x04, 0x92, 0x40, 0x04, 0xe2,
//		0x5c, 0x44, 0x92, 0x44, 0x38, 0xe2, 0x38
//	};
//	unsigned char *bitmaps[1] = { bitmap };
//	struct jbg_enc_state se;
//
//	jbg_enc_init(&se, 23, 5, 1, bitmaps,
//		output_bie, stdout);              /* initialize encoder */
//	jbg_enc_out(&se);                                    /* encode image */
//	jbg_enc_free(&se);                    /* release allocated resources */
//	printf("\n---------- end ----------\n");
//	system("pause");
//	return 0;
//}