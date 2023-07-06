#include <utility>
#include <set>
#include <cmath>

/*******************************************************************************/
/*** CREATION FUNCTIONS ***/

  void  readCells();
  void  readCollagen();
  void  putCells();
  void  putCellsHelper(int, int, int);
  void  calculatePerimeter(int);
  void  putCollagen();
  void  putCollagenHelper(int, int, double);

  double minDistance;
  double maxDistance;

/*******************************************************************************/
/*** Reads us in some cells ***/

void readCells()
{
  int x,y,t;
  FILE* inFile=fopen("lattice_.txt","r");
  if(inFile==NULL){
    printf("\n Nope... missing lattice_.txt\n\n");
    exit(0);
  }
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      fscanf(inFile,"%d",&lattice[i][j][0]);
      if(lattice[i][j][0]>numCells)
        lattice[i][j][0]-=COLLAGEN_OFFSET+PERIMETER_OFFSET*numCells;
      if(lattice[i][j][0]>0){
        cellVolumeList[lattice[i][j][0]].insert( std::make_pair(i, j) );
      }
    }
  }
  fclose(inFile);
  for(int cell=1;cell<=numCells; cell++)
    calculatePerimeter(cell);
}

/*******************************************************************************/
/*** Reads us in some collagen ***/

void readCollagen()
{
  int x,y,t;
  FILE* inFile=fopen("collagen_.txt","r");
  if(inFile==NULL){
    printf("\n Nope... missing collagen_.txt\n\n");
    exit(0);
  }
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      fscanf(inFile,"%d",&lattice[i][j][1]);
  fclose(inFile);
}

/*******************************************************************************/
/*** Populates an empty upper lattice ***/

void putCells()
{
    if(numCells==1)
        putCellsHelper(N/2,N/2,1);
    else if (numCells==2) {
        putCellsHelper(N/2-cellRadius,N/2-cellRadius-2,1);
        putCellsHelper(N/2+cellRadius,N/2+cellRadius+2,2);
    }
    else if (numCells==3) {
        putCellsHelper(N/2-2*cellRadius,N/2,1);
        putCellsHelper(N/2,N/2,2);
        putCellsHelper(N/2+2*cellRadius,N/2,3);
    }
    else if(numCells==7) {
        putCellsHelper(N/2-cellRadius,N/2-2*cellRadius,1);
        putCellsHelper(N/2+cellRadius,N/2-2*cellRadius,2);
        putCellsHelper(N/2-2*cellRadius,N/2,3);
        putCellsHelper(N/2,N/2,4);
        putCellsHelper(N/2+2*cellRadius,N/2,5);
        putCellsHelper(N/2-cellRadius,N/2+2*cellRadius,6);
        putCellsHelper(N/2+cellRadius,N/2+2*cellRadius,7);
    }
    else if(numCells==9) {
        putCellsHelper(N/2-2*cellRadius,N/2-2*cellRadius,1);
        putCellsHelper(N/2-2*cellRadius,N/2,2);
        putCellsHelper(N/2-2*cellRadius,N/2+2*cellRadius,3);
        putCellsHelper(N/2,N/2-2*cellRadius,4);
        putCellsHelper(N/2,N/2,5);
        putCellsHelper(N/2,N/2+2*cellRadius,6);
        putCellsHelper(N/2+2*cellRadius,N/2-2*cellRadius,7);
        putCellsHelper(N/2+2*cellRadius,N/2,8);
        putCellsHelper(N/2+2*cellRadius,N/2+2*cellRadius,9);
    }
    else {
      minDistance = 2;
      maxDistance = 100;
      
      int xCoordinate[numCells];
      int yCoordinate[numCells];
      int cell = 1;
      while (cell<=numCells){
          int x = rand()%(N/2) + N/4;
          int y = rand()%(N/2) + N/4;
          
          bool tooClose = false;
          bool tooFar = true;
          
          if (cell == 1)
              tooFar = false;
          
          for (int count=0; count<(cell-1); count++){
              int distanceSq = pow(x-xCoordinate[count],2) + pow(y-yCoordinate[count],2);
              int distance = sqrt(distanceSq);
              if (distance < cellRadius * minDistance)
                  tooClose = true;
              if (distance < cellRadius * maxDistance)
                  tooFar = false;
          }
          
          if (tooClose == false && tooFar == false) {
              putCellsHelper(x, y, cell);
              xCoordinate[cell-1] = x;
              yCoordinate[cell-1] = y;
              cell++;
          }
      }
  }
    
  for(int cell=1;cell<=numCells; cell++)
      calculatePerimeter(cell);
}

/*******************************************************************************/
/*** Places a single cell into the upper lattice ***/

void putCellsHelper(int x0, int y0, int cell)
{
  for(int x=-(int)(cellSpawn+0.5); x<=(int)(cellSpawn+0.5); x++){
    int ymax=(int)(sqrt(cellSpawn*cellSpawn-(double)(x*x)+0.5));
    for(int y=-ymax; y<=ymax; y++){
      if(lattice[(N+x0+x)%N][(N+y0+y)%N][0]==0){
		lattice[(N+x0+x)%N][(N+y0+y)%N][0]=cell;
		cellVolumeList[ cell ].insert( std::make_pair( (N+x0+x)%N , (N+y0+y)%N ) ); 
      }
    }
  }
}
/*******************************************************************************/
/*** Calculates cell perimeters ***/

void calculatePerimeter(int cell)
{
  for( std::set< std::pair<int, int> >::const_iterator it = cellVolumeList[cell].begin(); it!= cellVolumeList[cell].end(); ++it){ 
	int i = it->first;
	int j = it->second;
    if( lattice[i][j][0]!=lattice[(i+1)%N][j][0] ||
        lattice[i][j][0]!=lattice[i][(j+1)%N][0] || 
        lattice[i][j][0]!=lattice[(N+i-1)%N][j][0] || 
        lattice[i][j][0]!=lattice[i][(N+j-1)%N][0] )
    {
	  cellPerimeterList[ cell ].insert( std::make_pair(i, j) );
    }
  }
}

/*******************************************************************************/
/*** Puts collagen into the lower lattice ***/

void putCollagen()
{
  for(int n=0; n<numCollagen; n++){
    int x = rand()%N;
    int y = rand()%N;
    double m = tan((double)rand()/((double)RAND_MAX*2.0/3.14159));
    putCollagenHelper(x,y,m);
  }
}

/*******************************************************************************/
/*** Puts a single line of collagen into the lower lattice ***/

void putCollagenHelper(int x0, int y0, double slope)
{
  int sx = 2*(rand()%2)-1;
  int sy = 2*(rand()%2)-1;
  for(int i=0;i<collagenWidth;i++){
    int x  = x0+i*sy;
    int y  = y0-i*sx;
    double error = slope;
    int j=0;
    while(j<N){
      lattice[x][y][1]=1;
      while(error>0.5){
        y=(N+y+sy)%N;
        lattice[x][y][1]=1;
        error-=1.0;
        if(slope>1.0)
          j++;
      }
      x=(N+x+sx)%N;
      error+=slope;
      if(slope<=1.0)
        j++;
    }
  }
}

/*******************************************************************************/
