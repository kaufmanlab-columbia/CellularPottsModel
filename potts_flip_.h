#include <utility>
#include <set>
#include <vector>
#include <map>

/*******************************************************************************/


/*** SPIN FLIP FUNCTIONS ***/

    int    flip();
    void   choose();
    void	 addVolume( int i, int j, int cell);
    void	 removeVolume( int i, int j, int cell);
    void   adjustPerimeters(int cell);
    bool   maintainsContiguity();
    bool    noHoles();
    std::map< std::pair<int, int> , int > calculateChunkSites(int, int);

    int    iSite;
    int    jSite;
    int    oldCell;
    int    newCell;

/*******************************************************************************/
/*** Flip a spin and accept or reject ***/

int flip()
{

  // Choose a spin
  // The site to be invaded is stored in iSite, jSite
  // newCell invades oldCell

  choose();

  std::map< std::pair<int, int> , int > chunk;
  chunk = calculateChunkSites(iSite, jSite);

  // Subtract out parts of old energy affected by the flip chunk sites
  
    double deltaEnergy = 0.0;
    //double deltaBlobularEnergy = 0.0;
    double deltaPerimeterEnergy = 0.0;
    double deltaInteractionEnergy = 0.0;
    double deltaVolumeEnergy = 0.0;
    double deltaInternalEnergy = 0.0;

    
    for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
        std::pair<int,int> point = it->first;
        
        // the sites within the flip chunk itself
        deltaInteractionEnergy -= outplaneEnergy( point.first, point.second );
        deltaInteractionEnergy -= inplaneEnergy( point.first, point.second );
        
        // the four neighbor directions
        // if the neighbor is on the border of the chunk site, also subtract out its energy
        if( chunk.find( std::make_pair( (point.first + 1)%N, point.second ) ) == chunk.end() ){
            deltaInteractionEnergy -= inplaneEnergy( (point.first + 1)%N, point.second );
            deltaInteractionEnergy -= outplaneEnergy( (point.first + 1)%N, point.second );
        }
        if( chunk.find( std::make_pair( (N + point.first - 1)%N, point.second ) ) == chunk.end() ){
            deltaInteractionEnergy -= inplaneEnergy( (N + point.first - 1)%N, point.second );
            deltaInteractionEnergy -= outplaneEnergy( (N + point.first - 1)%N, point.second );
        }
        if( chunk.find( std::make_pair( point.first, (point.second + 1)%N ) ) == chunk.end() ){
            deltaInteractionEnergy -= inplaneEnergy( point.first, (point.second + 1)%N );
            deltaInteractionEnergy -= outplaneEnergy( point.first, (point.second + 1)%N );
        }
        if( chunk.find( std::make_pair( point.first, ( N + point.second - 1)%N ) ) == chunk.end() ){
            deltaInteractionEnergy -= inplaneEnergy( point.first, ( N + point.second - 1)%N );
            deltaInteractionEnergy -= outplaneEnergy( point.first, ( N + point.second - 1)%N );
        }
    }
    
  if(oldCell!=0){
      deltaVolumeEnergy -= volumeEnergy(oldCell);
 //   deltaEnergy -= anisotropyEnergy(oldCell);
 //   deltaBlobularEnergy -= blobularEnergy(oldCell);
      deltaPerimeterEnergy -= perimeterEnergy(oldCell);
      deltaInternalEnergy -= internalEnergy(oldCell);
  }

  if(newCell!=0){
    deltaVolumeEnergy -= volumeEnergy(newCell);
 //   deltaEnergy -= anisotropyEnergy(newCell);
 //   deltaBlobularEnergy -= blobularEnergy(newCell);
      deltaPerimeterEnergy -= perimeterEnergy(newCell);
      deltaInternalEnergy -= internalEnergy(newCell);

  }
    
  // Flip it

  for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
	  std::pair<int,int> point = it->first;
	  lattice[ point.first ][ point.second ][0]=newCell;
	  addVolume( point.first, point.second, newCell );
	  removeVolume( point.first, point.second, oldCell );
  }
  adjustPerimeters( newCell );
  adjustPerimeters( oldCell );


  // Add in energy associated with flipped site
    
    for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
        std::pair<int,int> point = it->first;
        deltaInteractionEnergy += outplaneEnergy( point.first, point.second );
        deltaInteractionEnergy += inplaneEnergy( point.first, point.second );
        
        // the four neighbor directions
        if( chunk.find( std::make_pair( (point.first + 1)%N, point.second ) ) == chunk.end() )
        {
            deltaInteractionEnergy += inplaneEnergy( (point.first + 1)%N, point.second );
            deltaInteractionEnergy += outplaneEnergy( (point.first + 1)%N, point.second );
        }
        
        if( chunk.find( std::make_pair( (N + point.first - 1)%N, point.second ) ) == chunk.end() ){
            deltaInteractionEnergy += inplaneEnergy( (N + point.first - 1)%N, point.second );
            deltaInteractionEnergy += outplaneEnergy( (N + point.first - 1)%N, point.second );
        }
        
        if( chunk.find( std::make_pair( point.first, (point.second+1)%N ) ) == chunk.end() ){
            deltaInteractionEnergy += inplaneEnergy( point.first, (point.second+1)%N );
            deltaInteractionEnergy += outplaneEnergy( point.first, (point.second+1)%N );
        }
        
        if( chunk.find( std::make_pair( point.first, ( N + point.second - 1)%N ) ) == chunk.end() ){
            deltaInteractionEnergy += inplaneEnergy( point.first, ( N + point.second - 1)%N );
            deltaInteractionEnergy += outplaneEnergy( point.first, ( N + point.second - 1)%N  );
        }
    }
 
    
  if(oldCell!=0){
    deltaVolumeEnergy += volumeEnergy(oldCell);
   // deltaEnergy += anisotropyEnergy(oldCell);
   // deltaBlobularEnergy += blobularEnergy(oldCell);
    deltaPerimeterEnergy += perimeterEnergy(oldCell);
      deltaInternalEnergy += internalEnergy(oldCell);

  }

  if(newCell!=0){
    deltaVolumeEnergy += volumeEnergy(newCell);
  //  deltaEnergy += anisotropyEnergy(newCell);
  //  deltaBlobularEnergy += blobularEnergy(newCell);
    deltaPerimeterEnergy += perimeterEnergy(newCell);
      deltaInternalEnergy += internalEnergy(newCell);

  }

   // deltaEnergy += deltaBlobularEnergy;
    deltaEnergy += deltaPerimeterEnergy;
    deltaEnergy += deltaVolumeEnergy;
    deltaEnergy += deltaInteractionEnergy;
    deltaEnergy += deltaInternalEnergy;

  // Accept or reject the flip
    
  if( deltaEnergy < 0 ){
    totalEnergy+=deltaEnergy;
      //totalBlobularEnergy+=deltaBlobularEnergy;
      totalPerimeterEnergy+=deltaPerimeterEnergy;
      totalVolumeEnergy += deltaVolumeEnergy;
      totalInteractionEnergy += deltaInteractionEnergy;
      totalInternalEnergy += deltaInternalEnergy;
    return 1;
  }
    
  else if( exp(-1.0*beta*deltaEnergy) > ((double) rand() / RAND_MAX)){
    totalEnergy+=deltaEnergy;
        //totalBlobularEnergy+=deltaBlobularEnergy;
        totalPerimeterEnergy+=deltaPerimeterEnergy;
        totalVolumeEnergy += deltaVolumeEnergy;
        totalInteractionEnergy += deltaInteractionEnergy;
        totalInternalEnergy += deltaInternalEnergy;

    return 1;
  }
  else{
	  for(std::map< std::pair<int, int> , int >::iterator it = chunk.begin(); it != chunk.end(); ++it){
		  std::pair<int,int> point = it->first;
		  int originalCellType = it->second;
		  lattice[ point.first ][ point.second ][0]=originalCellType;
		  if( originalCellType == oldCell ){
			  removeVolume( point.first, point.second, newCell );
			  addVolume( point.first, point.second, oldCell );
		  }
		  // if originalCellType was already newCell, then nothing has changed
		  // that is, the invasion did not cause this i,j site to change type
	  }
	  adjustPerimeters( newCell );
	  adjustPerimeters( oldCell );
	  return 0;
  }

}

/*******************************************************************************/
/*** Chooses a spin to flip ***/

void choose()
{

  do{

    int which,thing;

    // Choose a cell
    oldCell = (rand()%numCells)+1;

    // Choose a perimeter site within that cell to either extend or surrender
    std::set< std::pair<int, int> >::const_iterator it(cellPerimeterList[oldCell].begin());
    int p = rand()%cellPerimeterList[oldCell].size();
    advance(it,p);
    iSite = it->first;
    jSite = it->second;


    do{
      which = rand()%2;
      thing = 2*(rand()%2)-1;
      if(which==0)
         newCell = lattice[(N+iSite+thing)%N][jSite][0];
      else
        newCell = lattice[iSite][(N+jSite+thing)%N][0];
    }while( oldCell == newCell );

	  // Choose to invade or surrender
	  // (that is, if you randomly generate an even number, then
	  // switch from surrendering to invading)

    if(rand()%2==0){
      int tmp = oldCell;
      oldCell = newCell;
      newCell = tmp;
      if(which==0)
        iSite = (N+iSite+thing)%N;
      else
        jSite = (N+jSite+thing)%N;
    }

  }while(!maintainsContiguity());

  return;
}

/*******************************************************************************/
/*** Checks if the invasion would cause a cell to break into multiple pieces ***/
/*** Returns FALSE if so (returns TRUE if the flip would maintainContiguity) ***/

bool maintainsContiguity() {
    
  /*  if (oldCell == 0)
        return true; //there is no need to maintain contiguity for the environment*/
    
    
  std::vector<int> borders;

  // If the o's in the diagram below are sites in the chunk to the flipped
  // then borders are the x sites

  /* Diagram (chunkSize = 2 here, 0 marks (iSite,jSite) )
  x x x x x x x
  x o o o o o x
  x o o o o o x
  x o o 0 o o x
  x o o o o o x
  x o o o o o x
  x x x x x x x
  */

  /* add the top row of x's
  X X X X X X x
  x o o o o o x
  x o o o o o x
  x o o 0 o o x
  x o o o o o x
  x o o o o o x
  x x x x x x x
  */
  for( int i = iSite - chunkSize - 1; i <= iSite + chunkSize; i++ )
	  borders.push_back( lattice[ (N + i)%N ][ (N + jSite - chunkSize - 1)%N ][0]);

  /* add the right-side column of x's
  x x x x x x X
  x o o o o o X
  x o o o o o X
  x o o 0 o o X
  x o o o o o X
  x o o o o o X
  x x x x x x x
  */
  for( int j = jSite - chunkSize - 1; j <= jSite + chunkSize; j++ )
	  borders.push_back( lattice[ (iSite + chunkSize + 1)%N ][ (N + j)%N ][0] );

  // the bottom row
  for( int i = iSite + chunkSize + 1; i >= iSite - chunkSize; i-- )
	  borders.push_back( lattice[ (N + i)%N ][ (jSite + chunkSize + 1)%N ][0] );

  // the left-side column
  for( int j = jSite + chunkSize + 1; j >= jSite - chunkSize; j-- )
	  borders.push_back( lattice[ (N + iSite - chunkSize - 1)%N ][ (N + j)%N ][0]);

  // NOTE: THE BORDER SITES MUST BE ADDED IN CONTINUOUS ORDER FOR THIS
  // ALGORITHM TO WORK.


  /*
  int borders[8] = { lattice[(N+iSite-1)%N][(N+jSite-1)%N][0],
                     lattice[(N+iSite-1)%N][jSite][0],
                     lattice[(N+iSite-1)%N][(jSite+1)%N][0],
                     lattice[iSite][(jSite+1)%N][0],
                     lattice[(iSite+1)%N][(jSite+1)%N][0],
                     lattice[(iSite+1)%N][jSite][0],
                     lattice[(iSite+1)%N][(N+jSite-1)%N][0],
                     lattice[iSite][(N+jSite-1)%N][0] };
  */

  /*
	Count how many neighboring spins are in the same cell.
	If there are none, then this is the last spin of that cell.
	We will not let it disappear.
  */

  int totalCellCount = 0;
  for(int n=0; n<borders.size(); n++)
    if(borders[n] == oldCell)
      totalCellCount++;   
 
  if(totalCellCount==0)
     return false;

  /*
	The borders array can never be full of oldcell sites, so by
	starting at a non-cell site, you are guaranteed being able
	to traverse the entire contiguous region without breaks, as
	you aren't starting in the middle of a cell region.

	So, first find the first non-cell site.
   */

  int index = 0;
  while(index < borders.size() && borders[index] == oldCell)
    index++;
     
  /*
	Then move along until the site just before the next oldcell spin.
  */

  while(index < borders.size()-1 && borders[index+1] != oldCell)
    index++;

  /*
	Starting at the first cell spin, go around until you have
	hit every cell spin.  If you hit a non-cell spin while
	doing this, flipping Site would break contiguity.
	We will not let this happen.
  */

  int inCellCount = 0;
  while (inCellCount < totalCellCount)
  {
    index=(index+1)%(borders.size());
    if(borders[index] != oldCell)
      return false;
    inCellCount++;
  }

  return true;
}



/*******************************************************************************/
/*** Find the square chunk of sites around the chosen flip site corresponding to chunkSize ***/

std::map< std::pair<int, int>, int > calculateChunkSites( int i, int j ){
	std::map< std::pair<int, int>, int > chunk;
	for( int x=i-chunkSize; x<=i+chunkSize; x++ ){
		for( int y=j-chunkSize; y<=j+chunkSize; y++){
			std::pair<int, int> p;
			p.first = (x+N)%N;
			p.second = (y+N)%N;

			chunk.insert( std::pair< std::pair<int,int> , int>(p, lattice[ p.first ][ p.second ][0]) );
		}
	}
	return chunk;
}

/*******************************************************************************/
/*** Adjusts volume after a spin flip ***/

void removeVolume(int i, int j, int cell)
{
  if(cell!=0){
	  std::set< std::pair<int,int> >::iterator it = cellVolumeList[cell].find( std::make_pair(i,j) );
	  if( it != cellVolumeList[cell].end() )
		cellVolumeList[cell].erase( it );
  }
  return;
}

void addVolume(int i, int j, int cell){
  if(cell!=0)
	  cellVolumeList[cell].insert( std::make_pair(i,j) );
  return;
}


/*******************************************************************************/
/*** Adjusts perimeters after a spin flip ***/

// MUST BE DONE AFTER ADJUSTING VOLUMES
void adjustPerimeters( int cell ){

	// shorter way to do this: see if any of the chunk sites need to be added to the perimeter
	// run through all old perimeter sites and if they no longer need to be part of the perimeter, erase them from cellPerimeterList

	if( cell != 0 ){
		cellPerimeterList[cell].clear();
		  for( std::set< std::pair<int, int> >::const_iterator it = cellVolumeList[cell].begin(); it!= cellVolumeList[cell].end(); ++it){
			int i = it->first;
			int j = it->second;
			if(lattice[i][j][0] != cell )
					printf("\nproblem, we have a cell site that thinks it's not in the cell: (%d, %d)", i, j);
		    if( lattice[i][j][0]!=lattice[(i+1)%N][j][0] ||
		        lattice[i][j][0]!=lattice[i][(j+1)%N][0] ||
		        lattice[i][j][0]!=lattice[(N+i-1)%N][j][0] ||
		        lattice[i][j][0]!=lattice[i][(N+j-1)%N][0] )
		    {
			  cellPerimeterList[ cell ].insert( std::make_pair(i, j) );
		    }
		  }
	}
}
