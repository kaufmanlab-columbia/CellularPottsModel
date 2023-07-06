/*******************************************************************************/

/*** CELL ANALYSIS FUNCTIONS ***/

  void    measureCells();
  double  measureAnisotropy(int);

/*******************************************************************************/
/*** Measure cell properties ***/

void measureCells()
{
    double cellAnisotropy[numCells+1];
    for(int cell=1;cell<=numCells;cell++)
        cellAnisotropy[cell]=measureAnisotropy(cell);

    for(int n=0;n<4;n++){
        avg[n]=0.0;
        dev[n]=0.0;
    }
    for(int cell=1;cell<=numCells;cell++){
        avg[0]+=cellVolumeList[cell].size();
        avg[1]+=cellPerimeterList[cell].size();
        avg[2]+=cellAnisotropy[cell];
        avg[3]+=static_cast<double>(cellPerimeterList[cell].size()*cellPerimeterList[cell].size())/static_cast<double>(cellVolumeList[cell].size());
    }
    for(int n=0;n<4;n++)
        avg[n]/=numCells;
  
    for(int cell=1;cell<=numCells;cell++){
        dev[0]+=(cellVolumeList[cell].size()-avg[0])*(cellVolumeList[cell].size()-avg[0]);
        dev[1]+=(cellPerimeterList[cell].size()-avg[1])*(cellPerimeterList[cell].size()-avg[1]);
        dev[2]+=(cellAnisotropy[cell]-avg[2])*(cellAnisotropy[cell]-avg[2]);
        dev[3]+=(static_cast<double>(cellPerimeterList[cell].size()*cellPerimeterList[cell].size())/static_cast<double>(cellVolumeList[cell].size()))*((static_cast<double>(cellPerimeterList[cell].size()*cellPerimeterList[cell].size())/static_cast<double>(cellVolumeList[cell].size()))-avg[3]);
    }
    for(int n=0;n<4;n++) {
        dev[n]/=(numCells);
        dev[n]=sqrt(dev[n]);
    }
    return;
}

/*******************************************************************************/
/*** Measure anisotropy ***/


double measureAnisotropy(int cell)
{

  int    i,j,xi,xf,yi,yf,dx,dy,s;
  int    dist,distmax,distind[4];
  double slope,error;

  /*
    Go through every pair of perimeter spins, find the pair
    which are furthest apart, and record their positions.
  */

  distmax=0;

  for(std::set< std::pair<int, int> >::const_iterator it1 = cellPerimeterList[cell].begin(); it1!=cellPerimeterList[cell].end(); ++it1){
	  xi = it1->first;
	  yi = it1->second;

	  for(std::set< std::pair<int, int> >::const_iterator it2 = cellPerimeterList[cell].begin(); it2!=cellPerimeterList[cell].end(); ++it2){
      xf = it2->first;
      yf = it2->second;
      dx=xf-xi-N*(int)floor((float)(xf-xi)/(float)N+0.499);
      dy=yf-yi-N*(int)floor((float)(yf-yi)/(float)N+0.499);
      dist=dx*dx+dy*dy;
      if(dist>distmax){
        distmax=dist;
        distind[0]=xi;
        distind[1]=xf;
        distind[2]=yi;
        distind[3]=yf;
      }


    }
  }

  xi=distind[0];
  xf=distind[1];
  yi=distind[2];
  yf=distind[3];

  /*
     Find the midpoint.
  */

  dx = (xf-xi)-N*(int)floor((float)(xf-xi)/(float)N+0.499);
  dy = (yf-yi)-N*(int)floor((float)(yf-yi)/(float)N+0.499);

  int xm=(N+xi+dx/2)%N;
  int ym=(N+yi+dy/2)%N;

  /*
    Calculate the slope of the line perpendicular to
    the line between our two perimeter spins.
  */

  if(dy!=0)
    slope = -1.0*(double)dx/(double)dy;
  else
    slope = 2.0*(double)N;

  if(slope>=0.0)
    s=1;
  else
    s=-1;

  /*
    Start at the midpoint, and go out along the perpendicular
    until you find something not in the cell.
  */

  xi = xm;
  yi = ym;
  error = fabs(slope);

  if(lattice[xi][yi][0]!=cell)
    goto donei;

  while(1){
    if(lattice[xi][yi][0]!=cell){
      xi=(N+xi-1)%N;
      goto donei;
    }
    while(error>0.5){
      yi=(N+yi+s)%N;
      error-=1.0;
      if(lattice[xi][yi][0]!=cell){
        yi=(N+yi-s)%N;
        goto donei;
      }
    }
    xi=(xi+1)%N;
    error+=fabs(slope);
  }

  donei:

  xf = xm;
  yf = ym;
  error = fabs(slope);

  if(lattice[xf][yf][0]!=cell)
    goto donef;

  while(1){
    if(lattice[xf][yf][0]!=cell){
      xf=(xf+1)%N;
      goto donef;
    }
    while(error>0.5){
      yf=(N+yf-s)%N;
      error-=1.0;
      if(lattice[xf][yf][0]!=cell){
        yf=(N+yf+s)%N;
        goto donef;
      }
    }
    xf=(N+xf-1)%N;
    error+=fabs(slope);
  }

  donef:

  // Find the length of the perpendicular line.

  dx = xf-xi-N*(int)floor((float)(xf-xi)/(float)N+0.5);
  dy = yf-yi-N*(int)floor((float)(yf-yi)/(float)N+0.5);

  // Return

  double ani;

  if(dx==0 && dy==0)
    ani = (double)distmax;
  else
    ani = sqrt((double)distmax/(double)(dx*dx+dy*dy));

  return ani;
}


