/****************************************************************************
 *SPLITP:  Split data file into separate files determined by grid cell indices
 *         and species
 *Version:
 *         Unknown date, created by Henry Zhang
 *         12/28/2011, documented by Henry Zhang
 *         01/05/2012, re-implemented by Henry Zhang to take care NZ irreularity
 *Copyright:  Robert Grant and Henry Zhang, University of Alberta
 ***************************************************************************/
#include<search.h>
#include<string.h>
#include<strings.h>
#include<stdlib.h>
#include<stdio.h>
/*
 * LEN: The maximal length of a file name plus the null char at the end.
 */                                                     
#define LEN 17

/*
 *nmax: The maximum length of buffer, header
 */
#define nmax 2048

/*
 *INDEX: User Defined macro to get the output file index.
 */
#define INDEX(x,y,p,nh,nv,nz) ((((x)-1)*(nv)*(nz))+(((y)-1)*(nz))+((p)-1))

/*
 *splitp_: function prototype. Must end with "_" to comply with fortran 
 *         function call convention
 */
void splitp_(int *NHW,int *NHE,int *NVN,int *NVS, int *NZ, char *infile);


/***************************************************************************
 *splitp_:            implementation of splitp_ in C
 *Input Parameters:
 *                NHW--Index of horizontal west grid cell
 *                NHE--Index of horizontal east grid cell
 *                NVN--Index of vertical north grid  cell
 *                NVS--Index of vertical south grid  cell
 *                NZ --Number of species
 *                infile--input file name 
 *************************************************************************/                
void splitp_(int *NHW,int *NHE,int *NVN,int *NVS, int *NZ, char *infile){


     /*
      *src: 
      *head:  buffer to hold header string
      *buf:   a temporary character buffer
      *line:  char pointer for a line
      *prefix:file name prefix
      *ptr:   char pointer 
      *
      */
     char src[LEN],head[nmax],buf[nmax],*lines,*prefix,*ptr;

     /*
      *fin:   input file handle
      *fouts: array of output file handles
      */
     FILE *fin,*fouts[2000];

     /* position within the file to be read */
     fpos_t *pos;

     /*
      *i,j,k  : iteration indices for NH, NV, and NZ
      *len1   : temporary variable to hold the length of a line 
      *records: user defined multiplier to control how many chars to read once
      *items  : the actual number of chars read by fread
      *num    : total number of chars to be read once by fread
      *count  : temporary variable
      *cycles : =nh*nv*nz (horizontal grids x vertical grids x number of species
      */
     int i,j,k,len1,records=100,items=0,num=0,count=0,cycles;

     /* 
      *nx:     grid cell index for WE
      *ny:     grid cell index for NS
      *np:     current index for species
      *index1: index for output file units
      */
     int nx,ny,np,index1;

     /*
      *nh: number of horizontal grid cells
      *nv: number of vertial grid cells
      */
     int nh,nv;

     nh=(*NHE)-(*NHW)+1;
     nv=(*NVS)-(*NVN)+1;
     cycles=nh*nv*(*NZ);
     if(cycles>=2000){
        printf("error: number of grid cells exceed 400\n"); 
	return;
     }
     pos=(fpos_t*)calloc(1,sizeof(fpos_t));

     /*trim spaces at the end of character string*/
     for(i=0;i<LEN-1 && infile[i]!=' ';i++){
	 src[i]=infile[i];
     }
     if(infile[i]==' ')src[i]='\0';
     fin=fopen(src,"r");
     prefix=(char *)calloc((5+strlen(src)),sizeof(char));

     /*read the header line*/
     fgets(head,nmax,fin);

     /*
      *iterate through each grid cell for each single specie, and open
      *an output file handle for it and write the head line to it as well
      */
     count=0;
     for(i=(*NHW);i<=(*NHE);i++){
        for(j=(*NVN);j<=(*NVS);j++){
	   for(k=0;k<(*NZ);k++){
		snprintf(prefix,5+strlen(src),"%02d%02d%d%s",(i),(j),(k+1),src+1);
		fouts[count]=fopen(prefix,"w");
		if(fouts[count]==NULL)return;
		fprintf(fouts[count],"%s",head);
		count++;
	   }
     	}
     }

     /*Get the current position within the file*/
     fgetpos(fin,pos);

     /*Read a line into the buf*/
     fgets(buf,nmax,fin);

     /*Set the position back to pos*/
     fsetpos(fin,pos);

     /*Calculate the length of the line just read*/
     len1=strlen(buf);

     /*Total number of characters to be read once*/
     num=len1*cycles*records;
     lines=(char *)calloc(num,sizeof(char));
     if(lines==NULL){
	printf("error: failed to allocate memory\n");
	return ;
     }

     /*Read the whole file buffer by buffer, and then split it line by line*/ 
     items=0;
     do{
	items=fread(lines,sizeof(char),num,fin);
        /*Read one full buffer*/
        if(items==num){
                ptr=lines;
                count=0;
		while(count<items){
                      nx=10*(ptr[0]-'0')+ptr[1]-'0'-(*NHW)+1;
                      ny=10*(ptr[2]-'0')+ptr[3]-'0'-(*NVN)+1;
                      np=ptr[4]-'0';
		      index1=INDEX(nx,ny,np,nh,nv,*NZ);
                      if(index1>=0 && index1<cycles)
	              fwrite(ptr+16,sizeof(char),len1-16,fouts[index1]);
                   count+=len1;
                   ptr+=len1;
                }
         }
        /*Read less than full. Normally, the last read from the file*/
	else if(items>0 && items<num){
                ptr=lines;
                for(count=0;count<items;){
                        nx=10*(ptr[0]-'0')+ptr[1]-'0'-(*NHW)+1;
                        ny=10*(ptr[2]-'0')+ptr[3]-'0'-(*NVN)+1;
                        np=ptr[4]-'0';
			index1=INDEX(nx,ny,np,nh,nv,*NZ);
			if(index1>=0 && index1<cycles)
				fwrite(ptr+16,sizeof(char),len1-16,fouts[index1]);
			ptr=ptr+len1;
			count+=len1;
                }
	}
     }while(items!=0 && !feof(fin));

     /*Release resources to avoid memory leakage*/
     free(lines);
     free(prefix);
     free(pos);

     /*Close all opened output file units*/
     for(i=0;i<cycles;i++){
         fclose(fouts[i]);
    }

    /*Close the input file*/
    fclose(fin);
}

