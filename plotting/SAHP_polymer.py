import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import random 
import argparse


class Blob:
  """
  Stores the geometric and indexing properties of one polymer blob.
  """

  def __init__(self):
    """
    Initializes a blob at the origin with zero-valued properties.

    Arguments
        None

    Returns
        None
    """
    self.position = np.array([0.,0.,0.])
    self.rg = 0.
    self.ree = 0.
    self.i = 0
   
class Chain:
    """
    Represents a polymer chain composed of connected blobs.
    """

    def __init__(self):
        """
        Initializes an empty polymer chain at 300 K.

        Arguments
            None

        Returns
            None
        """
        self.blobs = []
        self.nblobs = 0
        self.temperature = 300.
                  
    def add_blob(self):
        """
        Adds one default-initialized blob to the chain.

        Arguments
            None

        Returns
            None
        """
        self.blobs.append(Blob())
        self.nblobs +=1
            
    def build_polymer(self, rg, ree, nblobs):
        """
        Builds a linear polymer from radii of gyration and end-to-end distances of blobs.

        Arguments
            rg : array-like
                Radius of gyration for each blob, in nanometers.
            ree : array-like
                End-to-end distance for each blob, in nanometers.
            nblobs : int
                Number of blobs in the polymer.

        Returns
            None
        """
        
        lastr = 0.
        for i in np.arange(nblobs):
          self.add_blob()
          self.blobs[i].rg = rg[i]
          self.blobs[i].ree = ree[i]
          if i>0:
            self.blobs[i].position[0] = lastr + (self.blobs[i-1].ree+self.blobs[i].ree)/2.0
          else:
            self.blobs[i].position[0] = 0.
          lastr = self.blobs[i].position[0]
          self.blobs[i].i = i
        self.contacts = np.zeros([nblobs,nblobs])
        self.contacts_distance = np.zeros([nblobs,nblobs])
  
    def check_ex_vol(self):
        """
        Checks whether any non-neighboring blobs overlap.

        Arguments
            None

        Returns
            overlap_detection : int
                1 if an excluded-volume overlap is detected; otherwise 0.
        """
      #whether the move is actually accepted depends only on whether there is an overlap 
        overlap_detection = 0
        for i in np.arange(self.nblobs):
          for j in np.arange(i+2, self.nblobs):
            separation_vector = self.blobs[i].position - self.blobs[j].position
            center_distance = np.sqrt(np.sum(separation_vector*separation_vector))
            if center_distance < self.blobs[i].rg + self.blobs[j].rg - 0.7*(self.blobs[i].rg + self.blobs[j].rg):
              overlap_detection = 1
              return overlap_detection
            
        return overlap_detection
    
    def check_bond_lengths(self):
        """
        Checks that neighboring blobs remain separated by their bond lengths.

        Arguments
            None

        Returns
            bond_length_error : int
                1 if an incorrect bond length is detected; otherwise 0.
        """
        # sanity check, confirms that every blob is separated from its neighbor by the correct amount
        
        bond_length_error = 0
        for i in np.arange(1,self.nblobs):
            j = i-1
            separation_vector = self.blobs[i].position - self.blobs[j].position
            center_distance = np.sqrt(np.sum(separation_vector*separation_vector))
            if np.abs(self.blobs[i].ree/2.0 + self.blobs[j].ree/2.0-center_distance)> 0.001:
              bond_length_error = 1
              print("WARNING:Bond length incorrect between ")
              # print(i,j,center_distance, self.blobs[i].ree + self.blobs[j].ree)
              return bond_length_error 
        return bond_length_error
        
    def sample(self, cutoff): 
        """
        Samples blob-blob contacts using an excess-distance cutoff.

        Arguments
            cutoff : float
                Maximum excess distance for defining a contact, in nanometers.

        Returns
            None
        """
        #Sample the contacts based on the cutoff

        for i in np.arange(self.nblobs):
            for j in np.arange(i+1, self.nblobs):

                separation_vector = self.blobs[i].position - self.blobs[j].position
                center_distance = np.sqrt(np.sum(separation_vector*separation_vector))

                self.contacts_distance[i][j] = center_distance - self.blobs[i].rg - self.blobs[j].rg
                #center_distance is the distance between centers of mass. Approximate "closest" distance between blobs by subtracting radius of gyration
                if center_distance - self.blobs[i].rg - self.blobs[j].rg < cutoff:
                    self.contacts[i][j]+=1

        
def mc_move(polymer,displacement_scale_nm):
  """
  Proposes and conditionally accepts a Monte Carlo move of the polymer.

  Arguments
      polymer : Chain
          Polymer chain whose coordinates will be updated.
      displacement_scale_nm : float
          Maximum displacement scale for the randomly proposed move.

  Returns
      None
  """
  #actually do the move
  #save pre move coordinates
  oldrs = np.zeros([polymer.nblobs,3]) 
  for j in np.arange(polymer.nblobs):
    oldrs[j] = polymer.blobs[j].position
  #randomly select blob to move
  blob_to_move = random.randint(1,polymer.nblobs-1) 
  if (blob_to_move > polymer.nblobs):
    print("Error: blob index too large")
  if (blob_to_move == 0):
    print("Error: blob index is 0")

  bondr = polymer.blobs[blob_to_move].position - polymer.blobs[blob_to_move -1].position #current bond vector to previous blob
  bondlength = polymer.blobs[blob_to_move].ree/2.0+ polymer.blobs[blob_to_move-1].ree/2.0 #bond length to previous blob should be half the sum of end to end distances
  
  #generate move 
  movex = displacement_scale_nm*(random.random()-0.5)
  movey = displacement_scale_nm*(random.random()-0.5)
  movez = displacement_scale_nm*(random.random()-0.5)
  mover = np.array([movex,movey,movez])

  #rescale to preserve bond length
  newbondr = bondr + mover
  newbondlength = np.sqrt(np.sum(newbondr**2))
  newbondr = newbondr*bondlength/newbondlength
  deltar = newbondr- bondr
  bondr = bondr + deltar
  
  #propagate move to rest of chain (move all subsequent blobs by the same amount)
  for j in np.arange(blob_to_move,polymer.nblobs):
    polymer.blobs[j].position = polymer.blobs[j].position + deltar
  if polymer.check_ex_vol() == 1:
    for j in np.arange(polymer.nblobs):
      polymer.blobs[j].position= oldrs[j]


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Simulates the SAHP and plots the blob-blob contacts.")

    parser.add_argument("--input_path", required=True, help="Path to the input data.") 
    parser.add_argument("--output_path", required=True, help="Path to the output data.")
    parser.add_argument("--nsteps", required=False, default=50000, help="Number of steps (50000).")
    parser.add_argument("--eqsteps", required=False, default=50000, help="Number of equilibration steps (50000).")
    parser.add_argument("--cutoff", required=False, default=0.55, help="Cutoff distance (0.55 nm).")
    parser.add_argument("--samplefreq", required=False, default=100, help="Sampling frequency (100).")
    parser.add_argument("--num_blobs", required=False, default=17, help="Number of blobs (17).")
    parser.add_argument("--displacement_scale_nm", required=False, default=0.5, help="Random displacement parameter in nanometers; each coordinate moves between -0.25 and +0.25 nm (0.5).")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    nsteps = int(args.nsteps)
    eqsteps = int(args.eqsteps)
    cutoff = float(args.cutoff)
    samplefreq = int(args.samplefreq)
    nblobs = int(args.num_blobs)
    displacement_scale_nm = float(args.displacement_scale_nm)

    polymer = Chain()

    # ree is first, rg is second in Angstroms of A66, Ree of N of the first resid of a blob and C of the second resid of a blob
    blob_parameters = {}
    blob_parameters['p1'] = [list(range(23,31)), 15.613676433677666, 6.547983042910457]
    blob_parameters['h1a'] = [list(range(31,39)), 15.2790805273261, 6.525860259946443]
    blob_parameters['s1'] = [list(range(39,42)), 8.01581052496582, 3.5757477375075237]
    blob_parameters['h1b'] = [list(range(42,48)), 11.396317632602312, 5.2006518670447095]
    blob_parameters['p2'] = [list(range(48,55)), 13.846447595665317, 6.178642713881808]
    blob_parameters['h2a'] = [list(range(55,64)), 14.59078265021237, 6.63168534602582]
    blob_parameters['s2'] = [[64,65], 3.3256334666297027, 3.6271323374954436]
    blob_parameters['h2b'] = [list(range(66,73)), 17.990940554253186, 6.038054381749419]
    blob_parameters['p3'] = [list(range(73,88)), 22.764660786255, 9.869292563833937]
    blob_parameters['h3a'] = [list(range(88,92)), 8.992071592283823, 4.54134437007044]
    blob_parameters['s3'] = [[92], 3.360232276711713, 1.7399349227334717]
    blob_parameters['h3b'] = [list(range(93,98)), 12.87225462725218, 5.53726881525724]
    blob_parameters['s4'] = [[98], 3.41693471871443, 1.737145991109951]
    blob_parameters['h3c'] = [list(range(99,104)), 13.151943529650659, 5.6396185296875405]
    blob_parameters['s5'] = [[104], 2.85548788435783, 1.7377817603329386]
    blob_parameters['h3d'] = [list(range(105,112)), 15.865919627573465, 6.164484046896214]
    blob_parameters['s6'] = [[112,113], 5.730063688414618, 4.015241597914009]    

    blob_list = ['p1', 'h1a', 's1', 'h1b', 'p2', 'h2a', 's2', 'h2b', 'p3', 'h3a', 's3', 'h3b', 's4', 'h3c', 's5', 'h3d', 's6']    
    
    sahp_output_filename="A66_parameritization_SAHP.txt" 
    
    ree = []
    rg = []

    # Divide by 10 to convert from Angstroms to nanometers
    for blob in blob_list:
      ree.append((blob_parameters[blob][1])/10)
      rg.append((blob_parameters[blob][2])/10)
  
    rs = np.zeros([nblobs,3])
    polymer.build_polymer(rg, ree, nblobs)

    #Equilibration
    for step in np.arange(eqsteps):
      mc_move(polymer,displacement_scale_nm)    

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4))
    
    #Production
    nsamples = 0
    for step in np.arange(nsteps):
      mc_move(polymer,displacement_scale_nm)
      if np.remainder(step,samplefreq)==0:
        polymer.sample(cutoff)
        nsamples+=1.0   

        #creates a histogram of contact distance
        hist, bins = np.histogram(polymer.contacts_distance, bins=14, range=(-4, 10))       

        center = ((bins[:-1] + bins[1:]) / 2)
        hist = [item * 1.0 / (nsamples*1.0) * 1.0 for item in hist]
        ax.plot(center, hist, lw=1)
        ax.axvline(x=cutoff,lw=1)   

    #Sanity check - make sure all the bond lengths are still constrained at the end of the simulation
    polymer.check_bond_lengths()    

    #saves the contact_length histrogram
    # figname = "%s/%sdist.pdf" % (output_path, "contact_coil_dist")
    # plt.savefig(figname, bbox_inches='tight')    
    
    #symmetrize contact matrix
    polymer.contacts = np.maximum( polymer.contacts, polymer.contacts.transpose())/nsamples
    polymer.contacts[polymer.contacts == 0] = 1  #the self intercation of each blob is considered 100%    

    ncol = 1
    nrow = 1
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4))
    vmax= 50
    
    contact_frequency = polymer.contacts*100
    contact_frequency_reduced = contact_frequency[:,[0,1,3,4,5,7,8,9,11,13,15]]
    contact_frequency_no_s_blobs = contact_frequency_reduced[[0,1,3,4,5,7,8,9,11,13,15],:]   

    np.savetxt("%s/%s" %(input_path, sahp_output_filename), contact_frequency_no_s_blobs)
    imgp = ax.imshow(contact_frequency_no_s_blobs, origin='lower', aspect='equal', cmap='Purples', vmin=0, vmax=vmax)

    box = ax.get_position()  # [x0, y0, width, height]
    cbar_width = 0.02
    cbar_pad = 0.01
    cbar_x = box.x1 + cbar_pad
    cbar_ax = fig.add_axes([cbar_x, box.y0, cbar_width, box.height])
    cbar = fig.colorbar(imgp, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label("Contact frequency (%)", fontsize=15, rotation=270, labelpad=20)

    ax.vlines(x=6-0.5, ymin=0-0.5, ymax=6-0.5, color="black", lw=4)
    ax.vlines(x=7-0.5, ymin=7-0.5, ymax=11-0.5, color="black", lw=4)
    ax.hlines(y=6-0.5, xmin=0-0.5, xmax=6-0.5, color="black", lw=4)
    ax.hlines(y=7-0.5, xmin=7-0.5, xmax=11-0.5, color="black", lw=4)
    ax.vlines(x=6-0.5, ymin=7-0.5, ymax=11-0.5, color="black", lw=4)
    ax.vlines(x=7-0.5, ymin=0-0.5, ymax=6-0.5, color="black", lw=4)
    ax.hlines(y=6-0.5, xmin=7-0.5, xmax=11-0.5, color="black", lw=4)
    ax.hlines(y=7-0.5, xmin=0-0.5, xmax=6-0.5, color="black", lw=4) 

    ax.xaxis.set_major_locator(plticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(plticker.MultipleLocator(1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ax.set_title('SAHP', loc='center', fontsize='20')

    figname = "%s/%s.pdf" % (output_path, sahp_output_filename[:-4])   

    plt.savefig(figname, bbox_inches='tight', dpi=100)
