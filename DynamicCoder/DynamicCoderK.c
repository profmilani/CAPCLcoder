// COMPILE: gcc DynamicCoderK.c ac.c ac.h -o Octree.exe -Ofast 
// EXECUTE: DynamicCoderK.exe <Nvx> <Volume> <Probabilities name> <decompositions> <frame>
// EXAMPLE: DynamicCoderK.exe 512 frame0001.bin prob_vet2.pr 9 1

#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "ac.h" 
#include <string.h> 
#include <time.h>

int main(int argc, char **argv)
{
    FILE* fp;
    
    /* Read inputs */
    if (argc != 6) 
	{
        printf("Usage: <executable> <Nvx> <Volume> <Probabilities name> <decompositions> <frame>.\n");
        return -1;
    } 
    
    int Nvx_orig = (int) atoi(argv[1]);
    unsigned char *vol = (unsigned char*) malloc(sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
    fp = fopen(argv[2],"rb");
    fread(vol, sizeof(unsigned char),Nvx_orig*Nvx_orig*Nvx_orig,fp);
    fclose(fp);
      
	int decompositions = (int) atoi(argv[4]);
	
	int frame = (int) atoi(argv[5]);
    
    /* Compute parameters */
    
    int Nvx_final = (int) pow(2,decompositions);
        
    int MAXdec = (int) (log(Nvx_orig)/log(2));

    /* Computer current transform */
	int *prob_vet2 = (int*) calloc(256*256,sizeof(int));

	int *freq = (int*) calloc(256,sizeof(int));
	
	unsigned char *map = (unsigned char*) malloc(sizeof(unsigned char)*256*4*4*16);
		
	int i, j, k;
	int max;
	int indMax = 0;
	unsigned char val;
	
	unsigned char *mapcurr; 
	int *prob_vet;

    clock_t tStartmp = clock();
	
	/* If the current frame is not the first one */
	if (frame>1)
	{		
		fp = fopen(argv[3],"rb");
		fread(prob_vet2, sizeof(int),256*256,fp);
		fclose(fp);
		
	    /* Compute transform on the previous probabilities */
		for (k=0; k<256; k++)
		{
			mapcurr = &(map[k*256]);
			prob_vet = &(prob_vet2[k*256]);
			
			mapcurr[0] = 0;

			val = 255;
			
			for (i=1; i<256; i++)
			{
				max = -1;
				for (j=1; j<256; j++)
				{
					if (prob_vet[j]>max)
					{
						max = prob_vet[j];
						indMax = j;
					}
				}
				prob_vet[indMax] = -1;
				mapcurr[indMax] = val--;
			}		
		}
		
		fp = fopen(argv[3],"rb");
		fread(prob_vet2, sizeof(int),256*256,fp);
		fclose(fp);
		
	    /* Compute frequencies on the previous probabilities */
		int sum = 0;
		for (i=1; i<256; i++)
			{
				for (j=0; j<256; j++)
				{
					freq[map[i+j*256]]+= prob_vet2[i+j*256];
					sum+= prob_vet2[i+j*256];
				}
			}		
		
		double rat = ((double)(16383-256))/((double) sum);
		
		for (i=0; i<256; i++)
		{
			freq[i] = (int) (((double)freq[i])*rat);

			if (freq[i]==0)
				freq[i] = 1;
		}	
				
	}
	else /* If the current frame is the first one */
	{	
		/* Initialize the transform as an identity */	
		for (i=0; i<256; i++)
		{
			for (j=0; j<256; j++)
			{
				map[i+j*256] = (unsigned char) i;
			}
		}
			
		/* Initialize frequencies as uniform */	
		for (i=0; i<256; i++)
		{
			freq[i] = 1;
		}
	}
	
	int time = (int)(clock() - tStartmp) *1000 / CLOCKS_PER_SEC;

    printf("Time initialization = %d\n", time);
	
	
    /* ENCODING PROCEDURE */
   
   /* Variables */
    int ndec;
    
    int x, y, z;
    
    unsigned char curr, out;
		
    unsigned char v0, v1, v2, v3, v4, v5, v6, v7;
    
    unsigned char *full_new, *temp;
    
    unsigned char *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7;

	/* Initialize arithmetic coder*/
    ac_encoder ace;
    ac_model acm;
											
	ac_encoder_init (&(ace),"foo");
	ac_model_init (&(acm), 256, freq, 1);
    
    /* Auxiliary volumes, variables and pointers*/
    unsigned char* FULL = (unsigned char*) malloc(sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
    memcpy(FULL,vol, sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
		
	unsigned char* FULLnew = (unsigned char*) malloc(sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
	
	int Nvx = Nvx_orig;
    int Nvxsq = Nvx*Nvx;
    int Nvx2 = Nvx/2;
    int Nvxsq2 = Nvxsq/4;
    	
    int *symbols = (int*) calloc(MAXdec,sizeof(int));
    
    unsigned char* s = (unsigned char*) malloc(sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig/8);
	unsigned char *ps = &(s[0]);

	unsigned char* PREC = (unsigned char*) calloc(Nvx_orig*Nvx_orig*Nvx_orig,sizeof(unsigned char));
    unsigned char *prec;
    prec = &(PREC[0]);

	int pos;
    
    clock_t tStart = clock();

	/* Loop on the levels of the Octree */
    for (ndec=0;ndec<MAXdec;ndec++)
    {	
		/* First block */
        c0 = &(FULL[0]);
        c1 = c0+1;
        c2 = c0+Nvx;
        c3 = c2+1;
        c4 = c0+Nvxsq;
        c5 = c4+1;
        c6 = c4+Nvx;
        c7 = c6+1;
        
        full_new = &(FULLnew[0]);
        
		/* Raster scan */
        for (z=0; z<Nvx; z+= 2)
        {
            for (y=0; y<Nvx; y+= 2)
            {
                for (x=0; x<Nvx; x+= 2)
                {
					/* Current block */
					curr = 128*(*c7)+64*(*c6)+32*(*c5)+16*(*c4)+8*(*c3)+4*(*c2)+2*(*c1)+(*c0);
				
					/* If the current block is non-empty*/
					if (curr >0)
					{
						/* Update charachteristics of the current block, useful for encoding future blocks */
						if (x<Nvx-2) 
							*(prec+1) += (2*(*c7)+(*c1));
										
						if (y<Nvx-2) 
							*(prec+Nvx2) += (8*(*c7)+4*(*c2));
						
						if (z<Nvx-2) 
							*(prec+Nvxsq2) += (128*(*c7)+64*(*c6)+32*(*c5)+16*(*c4));
						
						/* Transform current block */
						pos = curr+256*(*prec);
						*ps++= map[pos];
						
						/* Update probability */
						(*(prob_vet2+pos))++;
						
                        symbols[ndec]++;
						
						/* Update next level*/
						*full_new = (unsigned char)1;
					}
					else
					{
						/* Update next level*/
						*full_new = (unsigned char)0;
					}					
                   
                    full_new++;
                    prec++;

                    c0 += 2;
                    c1 += 2;
                    c2 += 2;
                    c3 += 2;
                    c4 += 2;
                    c5 += 2;
                    c6 += 2;
                    c7 += 2;
                   
                }
                c0 += Nvx;
                c1 += Nvx;
                c2 += Nvx;
                c3 += Nvx;
                c4 += Nvx;
                c5 += Nvx;
                c6 += Nvx;
                c7 += Nvx;
                
            }
            c0 += Nvxsq;
            c1 += Nvxsq;
            c2 += Nvxsq;
            c3 += Nvxsq;
            c4 += Nvxsq;
            c5 += Nvxsq;
            c6 += Nvxsq;
            c7 += Nvxsq;
                  
        }
        
		
		temp = FULL;
        FULL = FULLnew;
        FULLnew = temp;
		
		
        Nvx = Nvx/2;
        Nvxsq = Nvxsq/4;
        Nvx2 = Nvx2/2;
		Nvxsq2 = Nvxsq2/4;
    }
    
	/* Arithmetic encoding, following the inverse order of levels */
    int *cum = (int*) malloc(sizeof(int)*MAXdec);
    cum[0] = 0;
    
    for (ndec=1;ndec<MAXdec;ndec++)
        cum[ndec] = cum[ndec-1]+symbols[ndec-1];
    
    int shift;
    
    for (ndec=MAXdec-1;ndec>= 0;ndec--)
    {
        shift = cum[ndec];
        ps = &(s[0])+shift; 
		
        for (i=0;i<symbols[ndec];i++)
        {          
			out = *ps++;
            ac_encode_symbol (&(ace),&(acm), (unsigned char) out);
        }
    }
	
	/* Compute outputs*/	
    int enc_time = (int)(clock() - tStart) *1000 / CLOCKS_PER_SEC;
    printf("Time encoding = %d \n", enc_time);

	int TOTALE = ac_encoder_bits (&ace);
    printf("TOTALE = %d \n",TOTALE);
    
	
	ac_encoder_done (&(ace));
	ac_model_done (&(acm));
    
    
    
    /* DECODING PROCEDURE */
        
	/* Compute inverse transform */
    int c, l, f, d;
    unsigned char *inv_map = (unsigned char*) malloc(sizeof(unsigned char)*256*4*4*16);
    for (c=0; c<256; c++)
    {
        for (l=0; l<4; l++)
        {
            for (f=0; f<4; f++)
            {
                for (d=0; d<16; d++)
                {
                    inv_map[(int) map[c+256*l+1024*f+4096*d]+256*l+1024*f+4096*d] = c;
                }
            }
        }
    }
    
   	/* Initialize coders */
    ac_decoder acd;
	ac_decoder_init (&(acd),"foo");
	ac_model_init (&(acm), 256, freq, 1);
   
	/* Auxiliary volumes, variables and pointers*/
    unsigned char *full;
    memset(PREC, 0, sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
    memset(FULL, 0, sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
    memset(FULLnew, 0, sizeof(unsigned char)*Nvx_orig*Nvx_orig*Nvx_orig);
	FULL[0] = 1;
	prec = &(PREC[0]);

    Nvx = 2;
	Nvxsq = Nvx*Nvx;
	Nvx2 = Nvx/2;
	Nvxsq2 = Nvx2*Nvx2;
	
	clock_t tStart2 = clock();

	/* Loop on the levels of the Octree */
    for (ndec=decompositions-1;ndec>= 0;ndec--)
    {
		/* First block */
        c0 = &(FULLnew[0]);
        c1 = c0+1;
        c2 = c0+Nvx;
        c3 = c2+1;
        c4 = c0+Nvxsq;
        c5 = c4+1;
        c6 = c4+Nvx;
        c7 = c6+1;
        
        full = &(FULL[0]);
        
        /* Raster scan */
        for (z=0; z<Nvx; z+= 2)
        {            
            for (y=0; y<Nvx; y+= 2)
            {
                for (x=0; x<Nvx; x+= 2)
                {      
					/* If the current block is non-empty*/
                    if (*full>0)
                    {	    
						/* Decode current block */
						curr = ac_decode_symbol (&acd,&acm);
                       				
						/* Inverse transform current block */
						out = inv_map[curr+256*(*prec)];
                                        
						*c0 = (unsigned char)((out&1)>0);
						*c1 = (unsigned char)((out&2)>0);
						*c2 = (unsigned char)((out&4)>0);
						*c3 = (unsigned char)((out&8)>0);
						*c4 = (unsigned char)((out&16)>0);
						*c5 = (unsigned char)((out&32)>0);
						*c6 = (unsigned char)((out&64)>0);
						*c7 = (unsigned char)((out&128)>0);
                              
						/* Update charachteristics of the current block, useful for decoding future blocks */
						if (x<Nvx-2) 
							*(prec+1) += (2*(*c7)+(*c1));
										
						if (y<Nvx-2) 
							*(prec+Nvx2) += (8*(*c7)+4*(*c2));
						
						if (z<Nvx-2) 
							*(prec+Nvxsq2) += (128*(*c7)+64*(*c6)+32*(*c5)+16*(*c4));
					
                    }
                    else
					{
						*c0 = 0;
						*c1 = 0;
						*c2 = 0;
						*c3 = 0;
						*c4 = 0;
						*c5 = 0;
						*c6 = 0;
						*c7 = 0;
					}
                    
                    full++;
					prec++;

                    c0 += 2;
                    c1 += 2;
                    c2 += 2;
                    c3 += 2;
                    c4 += 2;
                    c5 += 2;
                    c6 += 2;
                    c7 += 2;
                }
                c0 += Nvx;
                c1 += Nvx;
                c2 += Nvx;
                c3 += Nvx;
                c4 += Nvx;
                c5 += Nvx;
                c6 += Nvx;
                c7 += Nvx;
            }
            c0 += Nvxsq;
            c1 += Nvxsq;
            c2 += Nvxsq;
            c3 += Nvxsq;
            c4 += Nvxsq;
            c5 += Nvxsq;
            c6 += Nvxsq;
            c7 += Nvxsq;
        }
        
        temp = FULL;
        FULL = FULLnew;
        FULLnew = temp;
        
        Nvx = Nvx*2;
        Nvxsq = Nvx*Nvx;
		Nvx2 = Nvx/2;
		Nvxsq2 = Nvx2*Nvx2;
    }
    
    
	
    /* Compute outputs*/	
	int dec_time = (int)(clock() - tStart2) *1000 / CLOCKS_PER_SEC;
    printf("Time decoding = %d \n", dec_time);
    printf("Time total = %d \n", enc_time+dec_time);
    
	ac_decoder_done (&(acd));
	ac_model_done (&(acm));
       
    int count = 0;
    for (i=0;i<Nvx_final*Nvx_final*Nvx_final;i++)
        count+= (int) FULL[i];
    printf("check: %d; points: %d",(int) memcmp(FULL,vol,Nvx_orig*Nvx_orig*Nvx_orig),count);
      
	  
	// fp = fopen("decvol.bin","wb");
	// fwrite(&(FULL[0]), sizeof(unsigned char), Nvx_final*Nvx_final*Nvx_final, fp);
	// fclose(fp);
		
	fp = fopen("prob_vet2.pr","wb");
	fwrite(&(prob_vet2[0]), sizeof(int), 256*256, fp);
	fclose(fp);
	  
	
    
}


// gcc CODIFICA11Kvid.c ac.c ac.h -o CODIFICA11Kvid.exe -Ofast
// CODIFICA11Kvid.exe 