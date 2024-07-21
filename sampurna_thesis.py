#### Module written for data generation and analysis for MPhil thesis: Sampurna Mukherjee
#Created 14.22 pm, 25-4-24
###imports##
import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#%matplotlib inline
import nglview as nv
import prolif as plf
## specific function imports:
#distance
from MDAnalysis.analysis import contacts
#creating dataframe from indices and parsing lists into dictionary
from collections import Counter
## for RSA calculation using DSSP
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from pathlib import Path # get pdb file name without extension
##################################Distance and Indices using MDAnalysis#####################################################
class distance_and_indices(object):
    def __init__(self,pdb_file_name): 
        self.u = mda.Universe(pdb_file_name)
        self.pdb_file_name = pdb_file_name
        self.pdb_name_without_ext = Path(pdb_file_name).stem
        #for hla_ab; we only want the heavy atoms
        self.hla_selection_heavy = ("(resid 1:276 or resid 462:561 or resid 1261:1269)and not name H*")
        self.ab_selection_heavy = ("(resid 632:853 or resid 952:1172)and not name H*")
        # for interactions with water, we want to keep hydrogens
        self.hla_selection = ("(resid 1:276 or resid 462:561 or resid 1261:1269)")
        self.ab_selection = ("(resid 632:853 or resid 952:1172)")
        self.water_selection =("name OW or name HW")
        self.whole_protein_selection = ("protein")
    def dist_calc_indice_generation(self,selection1,selection2):
        universe1 = self.u.select_atoms(selection1)
        universe2 = self.u.select_atoms(selection2)
        ##calculating distances
        dist_array = contacts.distance_array(universe1.positions,universe2.positions )
        #getting the indices of interface atoms
        self.interface_indices = np.where(dist_array<=4.5)
        return self.interface_indices
    def indice_to_df_static(self,selection1,selection2,ligand_name,interface_indices
                            ):
       #parsing the indices to create dataframe
        universe1 = self.u.select_atoms(selection1)
        universe2 = self.u.select_atoms(selection2)
        selection1_indices = interface_indices[0] #atom indices
        selection2_indices = interface_indices[1]
        # Getting resname and id from atom indices 
        selection1_list=[]
        selection2_list=[]
        selection1_atoml =[]
        selection2_atoml =[]
        resnames_list =[]
        for protein_atom, ab_atom in zip(selection1_indices,selection2_indices): 
            #The original function was for protein and ab, now it is made selection1 and 2 to include water.
            #I have changed some selections to selection1 and 2 instead og=f hla and ab, and same for universes
            #Too lazy to change the rest. SampurnaM, 08.02 am, 09-04-24
            
            #residue names
            selection1_res_name = universe1.atoms[protein_atom].resname #using the atom indice to get the corresponding resname
            selection2_res_name =universe2.atoms[ab_atom].resname
            selection1_atom_name =  universe1.atoms[protein_atom].name #.atoms.
            selection2_atom_name = universe2.atoms[ab_atom].name #.atoms
            #residue ids
            selection1_res_id = universe1.atoms[protein_atom].resid
            selection2_res_id = universe2.atoms[ab_atom].resid
            selection1_atom_id = protein_atom
            selection2_atom_id = ab_atom
        
            #combining both: resnameresid, like ASP90
            selection1_resnameid = selection1_res_name+str(selection1_res_id)
            selection2_resnameid = selection2_res_name+str(selection2_res_id)
            selection1_atomnameid = selection1_atom_name+str(selection1_atom_id)
            selection2_atomnameid = selection2_atom_name+str(selection2_atom_id)
            #appending this to a list
            selection1_list.append(selection1_resnameid)
            selection2_list.append(selection2_resnameid)
            selection1_atoml.append(selection1_atomnameid)
            selection2_atoml.append(selection2_atomnameid)
            #This will get a list of resname and resids of atoms. There will be (and there is) resid and resname repeaats
            #getting the counts of contacts per frame
        selection1_cnt_count = len(selection1_list)
        selection2_cnt_count = len(selection2_list)
        selection1_atom_count = len(selection1_atoml)
        selection2_atom_count = len(selection2_atoml)
        #getting the unique values and the count of times they have repeated using Counter
        #hla #using Counter
        unique_selection1_res = Counter(selection1_list).keys() 
        repeats_selection1_res = Counter(selection1_list).values() #gives the number of times each unique residue is in contact with another residue in the ab
        unique_selection1_at = Counter(selection1_atoml).keys() 
        repeats_selection1_at = Counter(selection1_atoml).values()
        #ab
        unique_selection2_res = Counter(selection2_list).keys() 
        repeats_selection2_res = Counter(selection2_list).values() #gives the number of times each unique residue is in contact with another residue in the hla
        unique_selection2_at = Counter(selection2_atoml).keys() 
        repeats_selection2_at = Counter(selection2_atoml).values()
        ## getting count of overall how many unique residues are there
        selection1_cnt_count_unique = len(unique_selection1_res)
        selection2_cnt_count_unique = len(unique_selection2_res)
    
        #creating array from this
        #processed_array = np.array((unique_hla_res,unique_ab_res)) # problem, as the values are of different size
        #so, instead, we can just return the two lists in a nested dictionary
        joined_dict = {
        "unique_protein_interface_res":unique_selection1_res, #didn't use f-strings on this because hard coded functions for this column down the class
        f"how_many_contacts_with_{ligand_name}_res": repeats_selection1_res,
        "unique_protein_name_interface_atoms":unique_selection1_at,
        f"how_many_contacts_with_{ligand_name}_atoms": repeats_selection1_at,
        f"unique_{ligand_name}_interface_res":unique_selection2_res, 
        f"how_many_contacts_with_protein_res": repeats_selection2_res,
        f"unique_{ligand_name}_interface_atoms":unique_selection2_at,
        f"how_many_contacts_with_protein_atoms": repeats_selection2_at
        
        #'unique_ab_contacts_pf': unique_ab_res,
        #'hla_res_in_cnt_with_ab':hla_cnt_count, 
        #'ab_res_in_cnt_with_hla':ab_cnt_count,
        #'unique_hla_contact_count':hla_cnt_count_unique,
        #'unique_ab_contact_count':ab_cnt_count_unique
        
                    }
        #returning the array
        #return joined_dict
          # joined dict has arrays of unequal length so it complains when I try to make it into df
        #solution: 
        #https://www.statology.org/pandas-dataframe-from-dict-with-different-length/
        indice_df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in joined_dict.items()]))
        return indice_df
    def indice_generation_hla_ab_water(self):
        #first generating the indices
        self.hla_ab_dist_output = self.dist_calc_indice_generation(self.hla_selection_heavy,self.ab_selection_heavy)
        self.hla_water_dist_output = self.dist_calc_indice_generation(self.hla_selection,self.water_selection)
        self.ab_water_dist_output = self.dist_calc_indice_generation(self.ab_selection,self.water_selection )
        ## then applying the indice_to_df_static function to generate dataframes
        self.hla_ab_indices = self.indice_to_df_static(self.hla_selection_heavy,self.ab_selection_heavy,"ab",self.hla_ab_dist_output)
        self.hla_water_indices = self.indice_to_df_static(self.hla_selection,self.water_selection,"water",self.hla_water_dist_output)
        self.ab_water_indices = self.indice_to_df_static(self.ab_selection,self.water_selection,"water",self.ab_water_dist_output)
############################################ ProLIF classes############################################################################
class prolif_generation_and_parsing(object):
    def __init__(self,pdb_file_name,indice_object):
        self.u = mda.Universe(pdb_file_name)
    def prolif_fp_generation(self,pickle_name,csv_name):
         #hla-ab
        hla_selection_prolif = ("resid 1:276 or resid 462:561 or resid 1261:1269")
        ab_selection_prolif = ("resid 632:853 or resid 952:1172")
        #water bridge
        water_selection_prolif = ("name OW or name HW")
       #creating the prolif universes
        hla_a1101_pdb_plf = self.u.select_atoms(hla_selection_prolif)
        ab_a1101_pdb_plf = self.u.select_atoms(ab_selection_prolif)
        water_a1101_pdb_plf = self.u.select_atoms(water_selection_prolif)
       #creating the  Prolif molecules
        self.hla_mol = plf.Molecule.from_mda(hla_a1101_pdb_plf)
        self.ab_mol = plf.Molecule.from_mda(ab_a1101_pdb_plf)
        self.water_mol = plf.Molecule.from_mda(water_a1101_pdb_plf)
        ##creating the fingerprint object
        self.hla_ab_fp = plf.Fingerprint( parameters={ "HBAcceptor": {"distance": 3.9,"DHA_angle":(90, 180)},
                                   "HBDonor": {"distance": 3.9,"DHA_angle":(90, 180)},
                                 "Cationic": {"distance": 4.0},
                                  "Anionic": {"distance": 4.0},
                                  },
           interactions = ['Anionic',
                     'CationPi',
                 'Cationic',
                 'EdgeToFace',
                 'FaceToFace',
                 'HBAcceptor',
                 'HBDonor',
                'Hydrophobic',
                 #'MetalAcceptor',
                 #'MetalDonor',
                 'PiCation',
                 'PiStacking',
                'VdWContact'
                 #'XBAcceptor',
                 #'XBDonor', 
                         ] , count=True
                    )
        self.water_hla_fp = plf.Fingerprint( parameters={ "HBAcceptor": {"distance": 3.9,"DHA_angle":(90, 180)},
                                   "HBDonor": {"distance": 3.9,"DHA_angle":(90, 180)},
                                 "Cationic": {"distance": 4.0},
                                  "Anionic": {"distance": 4.0},
                                  
                                 },
           interactions = ['Anionic',
                     'CationPi',
                 'Cationic',
                 'EdgeToFace',
                 'FaceToFace',
                 'HBAcceptor',
                 'HBDonor',
                'Hydrophobic',
                 #'MetalAcceptor',
                 #'MetalDonor',
                 'PiCation',
                 'PiStacking',
                'VdWContact'
                 #'XBAcceptor',
                 #'XBDonor', 
                         ] , count=True
                                            
                    )
        self.water_ab_fp =plf.Fingerprint( parameters={ "HBAcceptor": {"distance": 3.9,"DHA_angle":(90, 180)},
                                   "HBDonor": {"distance": 3.9,"DHA_angle":(90, 180)},
                                 "Cationic": {"distance": 4.0},
                                  "Anionic": {"distance": 4.0},
                                  
                                 },
           interactions = ['Anionic',
                     'CationPi',
                 'Cationic',
                 'EdgeToFace',
                 'FaceToFace',
                 'HBAcceptor',
                 'HBDonor',
                'Hydrophobic',
                 #'MetalAcceptor',
                 #'MetalDonor',
                 'PiCation',
                 'PiStacking',
                'VdWContact'
                 #'XBAcceptor',
                 #'XBDonor', 
                         ] , count=True
                    )

        
        #running the fingerprint object:hla-ab
        self.hla_ab_fp.run_from_iterable([self.hla_mol], self.ab_mol) ## HLA has less residues than Ab, so it is the ligand on this case
        self.hla_ab_ifp = self.hla_ab_fp.ifp
        #saving to pickle: using pickle_name parameter
        self.hla_ab_fp.to_pickle(pickle_name)
        #loading fingerprint back from pickle and converting to dataframe
        self.hla_ab_fp = plf.Fingerprint.from_pickle(pickle_name)
       
        #print("gOtcha!",self.hla_ab_ifp)
    # to troubeshoot,AttributeError: 'contacts_sm_static' object has no attribute 'hla_ab_fp'
        self.hla_ab_prolif_i_df = self.hla_ab_fp.to_dataframe()
        #running for hla- water*
        self.water_hla_fp.run_from_iterable([self.water_mol], self.hla_mol) #water is ligand here
        self.water_hla_ifp = self.water_hla_fp.ifp
        water_pickle_name = "hla_water"+pickle_name    
        self.water_hla_fp.to_pickle(water_pickle_name)
        #loading fingerprint back
        self.water_hla_fp= plf.Fingerprint.from_pickle(water_pickle_name)
        self.water_hla_prolif_i_df = self.water_hla_fp.to_dataframe()
        #running for ab-water
        self.water_ab_fp.run_from_iterable([self.water_mol], self.ab_mol) ## HLA has less residues than Ab, so it is the ligand on this case
        self.water_ab_ifp = self.water_ab_fp.ifp
        water_pickle_name1 = "ab_water"+pickle_name 
        self.water_ab_fp.to_pickle(water_pickle_name1)
        #loading fingerprint back
        self.water_ab_fp= plf.Fingerprint.from_pickle(water_pickle_name1)
        self.water_ab_prolif_i_df = self.water_ab_fp.to_dataframe()
             ############## Parsing starts#####################
    
        ##### getting the atom nmes: code from  https://github.com/chemosim-lab/ProLIF/issues/180
    def prolif_atom(self,fp_ifp,protein_mol):#prolif_df
        data = []
        index = []
        for i, frame_ifp in fp_ifp.items():
            index.append(i)
            frame_data = {}
            
            for (ligres, protres), residue_ifp in frame_ifp.items():
                for int_name, metadata_tuple in residue_ifp.items():
                    for metadata in metadata_tuple:
                        
                        # extract atom name using atom index in metadata and protein mol 
                        atoms = ",".join([ protein_mol.GetAtomWithIdx(x).GetMonomerInfo().GetName() for x in metadata["parent_indices"]["protein"] ]
                                        )
                        frame_data[(str(ligres), str(protres), int_name, atoms)] = 1
            
            data.append(frame_data) 
        ##now creating the dataframe from data
        out_df_name = pd.DataFrame( data,
                                      index=pd.Index(index, name="Frame")
                                  ) 
        out_df_name.columns = pd.MultiIndex.from_tuples(
                                out_df_name.columns, 
                                names=["ligand", "protein", "interaction", "atoms"]
                                                    )
        out_df_name = out_df_name.fillna(0).astype(int) 
        return out_df_name
        ################ now using this function to parse intermediate(i_dfs) to final prolif dfs
    def add_plf_atoms(self):
        
            
        self.hla_ab_prolif_df = self.prolif_atom(self.hla_ab_ifp,self.ab_mol)
        self.water_hla_prolif_df = self.prolif_atom(self.water_hla_ifp,self.hla_mol)
        self.water_ab_prolif_df = self.prolif_atom(self.water_ab_ifp,self.ab_mol)
    #### Addin NCIS to index df
    #### now we add NCIS to the prolif DF: function prolif_nci_add()
    def prolif_nci_add(self,prolif_df,index_df):
        pi_int_type = ['EdgeToFace','FaceToFace','PiCation','PiStacking','CationPi']
        res_list =[]
        res_int_dict = {}
        nci_category = ""
        # parse multilevel index 
        for i, level in zip(pd.MultiIndex.from_frame(prolif_df).names, 
                    pd.MultiIndex.from_frame(prolif_df).levels):
            
            #print(i,level)
            int_type = str(i[2])
            res_string = str(i[0])
            res = res_string.split(".")[0]
            #print(int_type,res)
            res_list.append(res)
            #print(res_list)
            if res in res_int_dict:
                res_int_dict[res].append(int_type)
            else:
                res_int_dict[res] = [int_type]
    ######## Now we just need to count number of NCIs per residue, i.e unique elements


            for j in index_df.unique_protein_interface_res:
                
                print("This is j",j)
                for key in res_int_dict.keys():
                    
                    print("This is key,",key)
                    if key==j:
                        
                        #print("This is key",key)
                        #print((Counter(res_int_dict[key])))
                        #print( sum( ( Counter( res_int_dict[key] ).values() ) ) )
                        #print(i,sum(Counter(i).values()))
           
                        val =sum(  Counter( res_int_dict[key] ).values() )
                        #sum( Counter( res_int_dict[key] ) )previous code was breaking the values into tuple, I had a for i in values loop
                        ## now getting the types of NCIs to make categories
                        for k in Counter(res_int_dict[key]).keys():
                            nci_category = nci_category+k
                            print(nci_category)
                            #print(k)
                             # Convert the string values to integers and then sum them up
                            #nci_sum = sum(int(val) for val in res_int_dict[key])
                            #nci_sum = sum(Counter(res_int_dict[key]).values()) 
                            #print("Value is ", nci_sum)
                            index_df.loc[index_df.unique_protein_interface_res == key,"NCI"] = val
                            index_df.loc[index_df.unique_protein_interface_res == key,"NCI_type"] = nci_category
                            nci_category ="" ### resetting the counter
        
        return index_df
##############################DSSP Modules##########################################################
class dssp_gen_and_parse(object):
    def __init__(self,pdb_file_name,pdb_name_without_ext,index_df_with_nci):
        self.pdb_file_name = pdb_file_name
        self.pdb_name_without_ext = pdb_name_without_ext
        self.index_df_with_nci = index_df_with_nci
        
    ##################### Now adding RSA information using Bio.DSSP#########
    def dssp_generation(self):
            
       ##using DSSP and generating the dssp object
        p = PDBParser(QUIET=True) # to silence the water warnings
        structure = p.get_structure(id = self.pdb_name_without_ext,
                                            file = self.pdb_file_name)
        model = structure[0]
        self.dssp = DSSP(model,self.pdb_file_name,dssp='mkdssp')
            ####### now parsing dssp to add RSA info in dataframes
    def dssp_parsing(self):
                
        #initializing empty list for residue index (indices) and RSA values: sasa
        self.indices =[]
        self.sasa =[]
        for data in self.dssp:
            index = data[0]
            self.indices.append(index)
            rsa = data[3]
            print("This is the RSA from generated DSSP dictionary",rsa)
            self.sasa.append(rsa)
                    
            ### now adding this informating by comparing with resname+id in dataframe 
            # updating code with .loc method of indexing
            #SampurnaM, 28-5-24
           
        for x in self.index_df_with_nci["unique_protein_interface_res"]:
            for y in zip(self.indices,self.sasa):
                dssp_resid_str = str(y[0])
                #brealking apart the unique residue name in df to separate digits and joining it back again
                resnum_list = [char for char in x if char.isdigit()]
                resnum = ''.join(resnum_list)
                #checking for match:
                if dssp_resid_str in resnum:
                    print("This is the RSA getting appended", y[1])
                    self.index_df_with_nci.loc[self.index_df_with_nci.unique_protein_interface_res == x,"RSA"] = y[1]
        return self.index_df_with_nci
############# adding CDR information################################################################
class cdr_annotation(object):
    def __init__(self,df_with_nci_and_sasa):
        self.df = df_with_nci_and_sasa
        
    def cdr_addition(self):
        #heavy
        cdr_h1 = range(657,665)
        cdr_h2 = range(682,690)
        cdr_h3 = range(727,743)
        #light
        cdr_l1 = range(977,986)
        cdr_l2 = range(1003,1010)
        cdr_l3 = range(1048,1057)
        self.df["CDR"] = ""
        #heavy cdrs
        for i in cdr_h1:
            i = str(i)
            mask = self.df.unique_ab_interface_res.str.contains(i)
            self.df.loc[mask,"CDR"] = "cdr_h1"
        for i in cdr_h2:
            i = str(i)
            mask = self.df.unique_ab_interface_res.str.contains(i)
            self.df.loc[mask,"CDR"] = "cdr_h2"
            
        for i in cdr_h3:
            i = str(i)
            mask = self.df.unique_ab_interface_res.str.contains(i)
            self.df.loc[mask,"CDR"] = "cdr_h3"
        ##light chain
        for i in cdr_l1:
            i = str(i)
            mask = self.df.unique_ab_interface_res.str.contains(i)
            self.df.loc[mask,"CDR"] = "cdr_l1"
        for i in cdr_l2:
            i = str(i)
            mask = self.df.unique_ab_interface_res.str.contains(i)
            self.df.loc[mask,"CDR"] = "cdr_l2"
        for i in cdr_l3:
            i = str(i)
            mask = self.df.unique_ab_interface_res.str.contains(i)
            self.df.loc[mask,"CDR"] = "cdr_l3"
        return self.df
