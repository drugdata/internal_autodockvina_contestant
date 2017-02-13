#!/usr/bin/env python


__author__ = 'j5wagner@ucsd.edu'

from d3r.celppade.custom_protein_prep import ProteinPrep
import os

chimera_prep_text = '''
import chimera
import sys
opened = chimera.openModels.open(sys.argv[1])
mol = opened[0]
import DockPrep
DockPrep.prep([mol])
from WriteMol2 import writeMol2
with open(sys.argv[2],'wb') as of:
    writeMol2([mol], of)
'''
class chimera_protprep(ProteinPrep):
    """Abstract class defining methods for a custom docking solution
    for CELPP
    """
    ProteinPrep.OUTPUT_PROTEIN_SUFFIX = '.pdbqt'
    def receptor_scientific_prep(self, 
                                 protein_file, 
                                 prepared_protein_file, 
                                 targ_info_dict={}):
        """
        Protein 'scientific preparation' is the process of generating
        a dockable representation of the candidate protein from a
        single-chain PDB file.
        :param protein_file: PDB file containing candidate protein.  
        :param prepared_protein_file: The result of preparation should have this file name.  
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.  
        :returns: True if preparation was successful. False otherwise.
        """
        with open('chimeraPrep.py', 'wb') as of:
            of.write(chimera_prep_text)
        os.system('grep ATOM ' + protein_file + ' > stripped_protein.pdb')
        os.system('chimera --nogui --script "chimeraPrep.py stripped_protein.pdb prepared_protein.mol2" 1> chimeraPrep.stdout 2> chimeraPrep.stderr')
        os.system('. /usr/local/mgltools/bin/mglenv.sh; $MGL_ROOT/bin/pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r prepared_protein.mol2 1> prepare_receptor4.stdout 2> prepare_receptor4.stderr')
        os.system('cp prepared_protein.pdbqt ' + prepared_protein_file)
        return True
        '''
        #Clean out all hetatms from original pdb
        data = open(protein_file).readlines()
        data = [line for line in data if not line[:6]=='HETATM']
        with open('stripped_protein.pdb','wb') as of:
            of.write(''.join(data))

        with open('chimeraPrep.py','wb') as of:
            of.write(chimera_prep_text) 
        os.system('chimera --nogui --script "chimeraPrep.py ' +
                  'stripped_protein.pdb prepared_protein.mol2' +
                  '" >& chimeraProtPrep.out')
        os.system('. /usr/local/mgltools/bin/mglenv.sh; pythonsh $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r prepared_protein.mol2')
        os.system('cp prepared_protein.pdbqt ' + prepared_protein_file)

        return True
        '''
        #return super(chimera_protprep,
        #             self).receptor_scientific_prep(protein_file, 
        #                                            prepared_protein_file, 
        #                                            targ_info_dict=targ_info_dict)
    


    
if ("__main__") == (__name__):
    import logging
    import os
    import shutil
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-p", "--pdbdb", metavar = "PATH", help = "PDB DATABANK which we will dock into")
    parser.add_argument("-c", "--challengedata", metavar="PATH", help = "PATH to the unpacked challenge data package")
    parser.add_argument("-o", "--prepdir", metavar = "PATH", help = "PATH to the output directory")
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level = logging.INFO )
    opt = parser.parse_args()
    pdb_location = opt.pdbdb
    challenge_data_path = opt.challengedata
    prep_result_path = opt.prepdir

    #running under this dir
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(prep_result_path), 'final.log')

    prot_prepper = chimera_protprep()
    prot_prepper.run_scientific_protein_prep(challenge_data_path, pdb_location, prep_result_path)

    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)
