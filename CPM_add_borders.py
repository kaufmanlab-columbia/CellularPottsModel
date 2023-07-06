import os
import numpy as np

def ID_perim_lattice(lattice,N=100):
    select_pix = []
    forbidden_pix = []
    # row iterator
    for nr,r in enumerate(lattice):
        # col iterator
        for nc,spin in enumerate(r):
            # up
            if spin!=lattice[nr][(nc+1)%N]:
                select_pix.append((nr,nc))
                forbidden_pix.append([spin,lattice[nr][(nc+1)%N]])

            # down
            elif spin!=lattice[nr][(nc-1)%N]:
                select_pix.append((nr,nc))
                forbidden_pix.append([spin,lattice[nr][(nc-1)%N]])

            # left
            elif spin!=lattice[(nr-1)%N][nc]:
                select_pix.append((nr,nc))
                forbidden_pix.append([spin,lattice[(nr-1)%N][nc]])

            # right
            elif spin!=lattice[(nr+1)%N][nc]:
                select_pix.append((nr,nc))
                forbidden_pix.append([spin,lattice[(nr+1)%N][nc]])
    
    
    return np.array(select_pix),np.array(forbidden_pix)

def second_pass(lattice,c,N=100):
    select_pix = []
    forbidden_pix = []
    
    # row iterator
    for nr,r in enumerate(lattice):
        # col iterator
        for nc,spin in enumerate(r):
            if spin==c+1:
                pass
            else:
            # up
                if spin!=lattice[nr][(nc+1)%N] and lattice[nr][(nc+1)%N]!=c+1:
                    select_pix.append((nr,nc))
                    forbidden_pix.append([spin,lattice[nr][(nc+1)%N]])

                # down
                elif spin!=lattice[nr][(nc-1)%N] and lattice[nr][(nc-1)%N]!=c+1:
                    select_pix.append((nr,nc))
                    forbidden_pix.append([spin,lattice[nr][(nc-1)%N]])

                # left
                elif spin!=lattice[(nr-1)%N][nc] and lattice[(nr-1)%N][nc]!=c+1:
                    select_pix.append((nr,nc))
                    forbidden_pix.append([spin,lattice[(nr-1)%N][nc]])

                # right
                elif spin!=lattice[(nr+1)%N][nc] and lattice[(nr+1)%N][nc]!=c+1:
                    select_pix.append((nr,nc))
                    forbidden_pix.append([spin,lattice[(nr+1)%N][nc]])
    
    
    return np.array(select_pix),np.array(forbidden_pix)

def perim_changer(lattice,p_info,cells=100):
    perim_x = []
    perim_y = []
    
    n_dict = {i:[] for i in range(0,cells+1)} # previously 1->cells+1
    
    for i in range(1,cells+1):
        first_1 = np.where(p_info[-1][:,0]==i)
        neighbors = np.unique(p_info[-1][first_1][:,1])
        for ii in neighbors:
            if ii in n_dict[i]:
                pass
            else:
                n_dict[ii].append(i)
                second_1 = np.where(p_info[-1][:,1]==ii)
                one = list(set(first_1[0]) & set(second_1[0]))
                x = list(p_info[0][one][:,0])
                y = list(p_info[0][one][:,1])
                perim_x.append(x)
                perim_y.append(y)

    perim_x_f = [j for i in perim_x for j in i]
    perim_y_f = [j for i in perim_y for j in i]
    
    lattice[tuple(np.vstack((perim_x_f,perim_y_f)))] = cells+1

    

def main(f_path,c,start_lat,end_lat,spacing=1000,N=100,fname='lattice_condensed.npy'):
    os.chdir(f_path)   
    all_lattice = np.load(fname)
    
    modified_lats = []
    for n in range(start_lat*N,(end_lat+1)*N,spacing*N):
        lat = all_lattice[n:n+N]
        p_info = ID_perim_lattice(lat,N=N)
        perim_changer(lat,p_info,cells=c)
        p_info = second_pass(lat,c,N=N)
        if p_info[0].shape[0]!=0:
            perim_changer(lat,p_info,cells=c)
        modified_lats.append(lat)
    
    modified_lats = np.vstack((modified_lats))
    np.save('lattice_select_perim.npy',modified_lats)
        



f_path = os.path.normpath(r"filepath")
os.chdir(f_path)
main(f_path,375,100,10000,spacing=100,N=150)
