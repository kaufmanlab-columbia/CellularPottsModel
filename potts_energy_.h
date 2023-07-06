/*******************************************************************************/
#include <set>
#include <utility>
#include <iostream>


/*** HAMILTONIAN FUNCTIONS ***/

    double  Hamiltonian();

    double  inplaneEnergy(int,int);
    double  outplaneEnergy(int,int);
    double  interactionEnergy(int);
    double  finalInteractionEnergy();

    double  volumeEnergy(int);
    double  finalVolumeEnergy();

    double  perimeterEnergy(int);
    double  finalPerimeterEnergy();

    double  siteInternalEnergy(int,int);
    double  internalEnergy(int);
    double  finalInternalEnergy();

    double  blobularEnergy(int);
    double  finalBlobularEnergy();
    double  anisotropyEnergy(int);
    double  measureAnisotropy( int );

/*******************************************************************************/
/*** Returns the in-plane interaction energy of a lattice site ****/

double inplaneEnergy(int a, int b)
{
  double energy = 0.0;

  // if the site is a cell site
  if( lattice[a][b][0] != 0 ){

	// if the site's upper neighbor is not the same type as itself
    if( lattice[a][b][0] != lattice[(a+1)%N][b][0] ){
      // if the site's upper neighbor is a cell
      if( lattice[(a+1)%N][b][0] > 0 )
        energy += J_cel;
      // if the site's upper neighbor is air
      else
        energy += J_air;
    }

    // lower neighbor
    if( lattice[a][b][0] != lattice[(N+a-1)%N][b][0] ){
      if( lattice[(N+a-1)%N][b][0] > 0 )
        energy += J_cel;
      else
        energy += J_air;
    }

    // right neighbor
    if( lattice[a][b][0] != lattice[a][(b+1)%N][0] ){
      if( lattice[a][(b+1)%N][0] > 0 )
        energy += J_cel;
      else
        energy += J_air;
    }

    // left neighbor
    if( lattice[a][b][0] != lattice[a][(N+b-1)%N][0] ){
      if( lattice[a][(N+b-1)%N][0] > 0 )
        energy += J_cel;
      else
        energy += J_air;
    }

  }

  return energy;
}

/*******************************************************************************/
/*** Returns the out-of-plane interaction energy of a lattice site ****/

double outplaneEnergy(int a, int b)
{
  if( lattice[a][b][0]!=0 && lattice[a][b][1]!=0 )
    return J_col;
  else
    return 0.0;
}

/*******************************************************************************/
/*** Returns the volume energy of a cell ***/

double volumeEnergy(int cell)
{
    return L_vol*((double)cellVolumeList[cell].size()-targetVolume)*((double)cellVolumeList[cell].size()-targetVolume);
}

/*******************************************************************************/
/****Return the perimeter energy *******/

double perimeterEnergy (int cell) {
    return L_per*((double)cellPerimeterList[cell].size() - targetPerimeter)*((double)cellPerimeterList[cell].size() - targetPerimeter);
}


/*******************************************************************************/
/*** Returns the interaction energy of a cell ****/

double interactionEnergy(int cell)
{
    double energy = 0.0;
    
    for(std::set< std::pair<int, int> >::iterator it = cellPerimeterList[cell].begin(); it!=cellPerimeterList[cell].end(); ++it)
        energy += inplaneEnergy( it->first , it-> second );
    
    for(std::set< std::pair<int, int> >::iterator it = cellVolumeList[cell].begin(); it!=cellVolumeList[cell].end(); ++it)
        energy += outplaneEnergy( it->first , it-> second );
    
    return energy;
}

/*******************************************************************************/
/*** Return the full perimeter energy of the lattice ***/

double finalPerimeterEnergy() {
    double energy = 0.0;
    for (int cell=1; cell<=numCells; cell++){
        energy += perimeterEnergy(cell);
    }
    return energy;
}

/*******************************************************************************/
/*** Return the full interaction energy of the lattice ***/

double finalInteractionEnergy() {
    double energy = 0.0;
    for (int cell=1; cell<=numCells; cell++){
        energy += interactionEnergy(cell);
    }
    return energy;
}

/*******************************************************************************/
/*** Return the full volume energy of the lattice ***/

double finalVolumeEnergy() {
    double energy = 0.0;
    for (int cell=1; cell<=numCells; cell++){
        energy += volumeEnergy(cell);
    }
    return energy;
}

/*******************************************************************************/
/*** Returns the full hamiltonian of the lattice ***/

double Hamiltonian()
{
  double energy = 0.0;
  for(int cell=1;cell<=numCells;cell++){
      energy += interactionEnergy(cell);
      energy += volumeEnergy(cell);
   //   energy += anisotropyEnergy(cell);
   //   energy += blobularEnergy(cell);
      energy += perimeterEnergy(cell);
      energy += internalEnergy(cell);
  }
  return energy;
}

/*******************************************************************************/
/*****return the internal energy of one site ********/


double siteInternalEnergy(int a, int b) {
    double energy = 0.0;
    
    //if the site is a cell site
    if(lattice[a][b][0] != 0){

        //if the site's upper neighbor is the same type as itself
        if(lattice[a][b][0] == lattice[(a+1)%N][b][0]){
            energy -= L_int;
        }
        //lower neighbor
        if(lattice[a][b][0] == lattice[(N+a-1)%N][b][0]){
            energy -= L_int;
        }
        //right neighbor
        if(lattice[a][b][0] == lattice[a][(b+1)%N][0]){
            energy -= L_int;
        }
        //left neighbor
        if(lattice[a][b][0] == lattice[a][(N+b-1)%N][0]){
            energy -= L_int;
        }

    }
    return energy;
}

/*******************************************************************************/
/*****return the internal energy of one cell********/


double internalEnergy(int cell) {
    double energy = 0.0;
    
    for(std::set<std::pair<int,int> >::iterator it = cellVolumeList[cell].begin(); it!=cellVolumeList[cell].end(); ++it){
        energy += siteInternalEnergy(it->first,it->second);
    }
    
    return energy;
}


/*******************************************************************************/
/*** Return the full internal energy of the lattice ***/

double finalInternalEnergy() {
    double energy = 0.0;
    for (int cell=1; cell<=numCells; cell++){
        energy += internalEnergy(cell);
    }
    return energy;
}

/*******************************************************************************/
/*** Return the full blobular energy of the lattice ***/


double finalBlobularEnergy() {
    double energy = 0.0;
    for (int cell=1; cell<=numCells; cell++){
        energy += blobularEnergy(cell);
    }
    return energy;
}

 
/*******************************************************************************/
/*** Returns the perimeter energy of a cell ***/

double anisotropyEnergy(int cell)
{
	//	return L_ani*(double)cellPerimeterList[cell].size()/(double)cellVolumeList[cell].size();
		return L_ani*measureAnisotropy( cell );
}

/*******************************************************************************/
/*** Returns the blobular energy of a cell ***/

double blobularEnergy(int cell)
{
    
    // (N is the total number of perimeter sites)
    // iterate over logN random perimeter site, and for each
    // iterate over logN random perimeter sites
    // and if the line between the site and the other site goes outside the cell
    // penalize
    
    int N = cellPerimeterList[cell].size();
    int logN = log( (double)N );
    int advanceAmount = logN;
    
    int energy = 0;
    int number = 0;
    
    int outerCount = 0;
    
    std::set< std::pair<int, int> >::const_iterator it1(cellPerimeterList[cell].begin());
    
    do{
        advance(it1, advanceAmount);
        outerCount += advanceAmount;
        
        int ai = it1->first;
        int aj = it1->second;
        
        int count = 0;
        std::set< std::pair<int, int> >::const_iterator it2(cellPerimeterList[cell].begin());
        do{
            advance(it2,advanceAmount);
            count += advanceAmount;
            
            int bi = it2->first;
            int bj = it2->second;
            
            int dx = bi-ai-N*(int)floor((float)(bi-ai)/(float)N+0.5);
            int sx = (dx>0)-(dx<0);
            int dy = bj-aj-N*(int)floor((float)(bj-aj)/(float)N+0.5);
            int sy = (dy>0)-(dy<0);
            
            double slope;
            if(dx!=0)
                slope = (double)dy/(double)dx;
            else
                slope = (double)N;
            
            int x = ai;
            int y = aj;
            double error = fabs(slope);
            
            do{
                number++;
                if(lattice[x][y][0]!=cell){
                    energy++;
                    goto done;
                }
                while(error>0.5){
                    y=(y+sy+N)%N;
                    number++;
                    if(lattice[x][y][0]!=cell){
                        energy++;
                        goto done;
                    }
                    error=error-1.0;
                }
                x=(x+sx+N)%N;
                error+=fabs(slope);
            }while(x!=(ai+dx+N)%N);
            
        done:;
        }while( count + advanceAmount < N );
        
    }while( outerCount + advanceAmount < N );
    
    
    return L_blb * (double)energy * log((double)energy) / (double)(cellPerimeterList[cell].size());
    
}

/*******************************************************************************/
