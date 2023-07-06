/*******************************************************************************/
/*** PRINTING FUNCTIONS ***/

    void  printLog(char*);
    void  printCells(char*);
    void  printLattice(char*);
    void  printFinalLattice(char*);
    void  printCollagen(char*);
    void  printRatioList(char*_);

/*******************************************************************************/
/*** Prints a log file ***/

void printLog(char *fname)
{
  FILE *pfile;
  pfile=fopen(fname,"w");
  fprintf(pfile,"SEED\t\t%d\n\n",seed);
  fprintf(pfile,"LOOPS\t\t%d\n",numLoops);
  fprintf(pfile,"FLIPS\t\t%d\n",numFlips);
  fprintf(pfile,"PRINT\t\t%d\n\n",numPrint);
  fprintf(pfile,"CHUNKSIZE\t\t%d\n\n",chunkSize);
  fprintf(pfile,"BETA\t\t%lf\n\n",beta);
  fprintf(pfile,"AIR\t\t%lf\n",J_air);
  fprintf(pfile,"CELL\t\t%lf\n\n",J_cel);
  fprintf(pfile,"COLLAGEN\t%lf\n\n",J_col);
  fprintf(pfile,"VOLUME\t\t%lf\n",L_vol);
  fprintf(pfile,"ANISOTROPY\t%lf\n",L_ani);
  fprintf(pfile,"INTERNAL\t%lf\n",L_int);
  fprintf(pfile,"PERIMETER\t%lf\n\n",L_per);
  fprintf(pfile,"LATTICE EDGE\t%d\n",N);
  fprintf(pfile,"NUM CELLS\t%d\n",numCells);
  fprintf(pfile,"NUM COLLAGEN\t%d\n",numCollagen);
  fprintf(pfile,"CELL SPAWN\t%lf\n",cellSpawn);
  fprintf(pfile,"CELL RADIUS\t%lf\n",cellRadius);
  fprintf(pfile,"COLLAGEN WIDTH\t%d\n\n",collagenWidth);
  fclose(pfile);
}

/*******************************************************************************/
/*** Prints a list Vol/Per ratio file ***/

void printRatioList(char *fname)
{
    FILE *pfile;
    pfile=fopen(fname,"w");
    fprintf(pfile,"%5s ","Number");
    fprintf(pfile,"%8.8s ","P^2/V");
    fprintf(pfile,"\n");
    
    for(int cell=1; cell<=numCells; cell++) {
        fprintf(pfile,"%5d ",cell);
        fprintf(pfile,"%8.3lf ", static_cast<double>(cellPerimeterList[cell].size()*cellPerimeterList[cell].size())/static_cast<double>(cellVolumeList[cell].size()));
        fprintf(pfile,"\n");
    }
    
    fclose(pfile);
}

/*******************************************************************************/

/*** Prints the cells to a text file ***/

void printCells(char *fname)
{
  FILE *printTo;
  printTo=fopen(fname,"w");
  for(int cell=1;cell<=numCells;cell++){
    fprintf(printTo,"%d\n", (int)cellVolumeList[cell].size());
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        if(lattice[i][j][0]==cell)
          fprintf(printTo,"%9d%9d%9d\n",i,j,2);
      }
    }
  }
  fclose(printTo);
}

/*******************************************************************************/
/*** Prints the lattice to a text file ***/

void printLattice(char *fname)
{
  FILE *pfile;
  pfile=fopen(fname,"w");
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      //air
      if(lattice[i][j][0]==0){
        if(lattice[i][j][1]==0)
          fprintf(pfile,"%d ",0);
        else
          fprintf(pfile,"%d ",numCells+COLLAGEN_OFFSET);
      }
      //cell
      else{
    	for(std::set< std::pair<int, int> >::const_iterator it = cellPerimeterList[ lattice[i][j][0] ].begin(); it!=cellPerimeterList[ lattice[i][j][0] ].end(); ++it){
    		if( (it->first == i) && (it->second == j) ){
            fprintf(pfile,"%d ",lattice[i][j][0]+COLLAGEN_OFFSET+PERIMETER_OFFSET*numCells); //add 1 to perimeter sites
            goto yep;
          }
        }
        fprintf(pfile,"%d ",lattice[i][j][0]);
        yep:;
      }
    }
    fprintf(pfile,"\n");
  }
  fclose(pfile);
}

/*******************************************************************************/
/*** Prints final lattice to a text file ***/

void printFinalLattice(char *fname)
{
    FILE *pfile;
    pfile=fopen(fname,"w");
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            //air
            if(lattice[i][j][0]==0){
                if(lattice[i][j][1]==0)
                    fprintf(pfile,"%d ",0);
                else
                    fprintf(pfile,"%d ",numCells+COLLAGEN_OFFSET);
            }
            //cell
            else
                fprintf(pfile,"%d ",lattice[i][j][0]);
        }
        fprintf(pfile,"\n");
    }
    fclose(pfile);
}

/*******************************************************************************/
/*** Prints the lower lattice only to a text file ***/

void printCollagen(char *fname)
{
  FILE *pfile;
  pfile=fopen(fname,"w");
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++)
      fprintf(pfile,"%d ",lattice[i][j][1]);
    fprintf(pfile,"\n");
  }
  fclose(pfile);
}

/*******************************************************************************/

