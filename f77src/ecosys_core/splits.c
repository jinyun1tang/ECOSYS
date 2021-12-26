#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define nmax 2048
void splits_(const int *NHW,const int *NHE,const int *NVN,const int *NVS, const char *outdir, const char *infile, int *failure)
{
     char srcl[256],head[nmax],buf[nmax],*lines,*prefix,*ptr;
     char src[17];
     char *outdirl;
     FILE *fin,*fouts[400];
     fpos_t *pos;
     int nh=(*NHE)-(*NHW)+1;
     int nv=(*NVS)-(*NVN)+1;
     int i,j,k,offset,len,records=100,items=0,num=0,count=0,cycles;
     static const char modfile[]=__FILE__;

     *failure=0;
     cycles=nh*nv;
     pos=(fpos_t*)malloc(sizeof(fpos_t));

     for(i=0;i<16 && infile[i]!=' ';i++){
   	   src[i]=infile[i];
     }
     if(infile[i]==' ')src[i]='\0';

     outdirl=(char *)malloc(sizeof(char)*(1+strlen(outdir)));

     for(i=0;i<strlen(outdir) && outdir[i]!=' ';i++){
   	   outdirl[i]=outdir[i];
     }
     outdirl[i]='\0';

     sprintf(srcl,"%s/%s", outdirl,src);
     fin=fopen(srcl,"r");
     if(fin==NULL){
       printf("Fail to open file %s in %s at line %d\n",srcl,modfile,__LINE__);

       *failure=1;
       return;
     }
     prefix=(char *)malloc(sizeof(char)*(5+strlen(src)));
     fgets(head,nmax,fin);
     k=0;
     /*for(k=0;k<cycles;k++){*/
     for(i=(*NHW);i<=(*NHE);i++){
      for(j=(*NVN);j<=(*NVS);j++){
/*	i=(k/(*nv));
	j=(k%(*nv));*/
	      snprintf(prefix,5+strlen(src),"%02d%02d%s",i,j,src);
        sprintf(srcl,"%s/%s",outdirl,prefix);
	      fouts[k]=fopen(srcl,"w");
	      if(fouts[k]==NULL){
          *failure=1;
          printf("Faile to create file %s in %s at %d\n",srcl, modfile, __LINE__);
          return;
        }
	      fprintf(fouts[k],"%s",head);
        k++;
     /*}*/
      }
     }
     fgetpos(fin,pos);
     fgets(buf,nmax,fin);
     fsetpos(fin,pos);
     len=strlen(buf);
     num=len*cycles*records;
     lines=(char *)malloc(sizeof(char)*num);
     if(lines==NULL){
       *failure=1;
	     printf("error: failed to allocate memory while splitting soil output in %s at line %d\n",modfile,__LINE__);
	     return;
     }
     items=0;
     offset=0;
     do{
	      items=fread(lines,sizeof(char),num,fin);
        if(items==num){
		      ptr=lines;
		      for(i=0;i<records;i++){
	          for(j=0;j<cycles;j++){
			        fwrite(ptr+16,sizeof(char),len-16,fouts[j]);
			        ptr=ptr+len;
			        offset++;
		        }
		      }
        }
      	else if(items>0 && items<num){
          ptr=lines;
		      count=0;
		      for(i=0;i<records;i++){
	          for(j=0;j<cycles;j++){
			        fwrite(ptr+16,sizeof(char),len-16,fouts[j]);
			        ptr=ptr+len;
			        count+=len;
			        offset++;
		        }
     		   if(count>=items)break;
		     }
	     }
     }while(items!=0 && !feof(fin));
     free(lines);
     free(prefix);
     free(pos);
     free(outdirl);
     for(i=0;i<cycles;i++){
         fclose(fouts[i]);
     }
    fclose(fin);
    return;
 }
