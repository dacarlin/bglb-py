import pyrosetta 
import rosetta 

mutant_name = 'H178K' 

with open( 'input_files/flags' ) as fn:
    flags = fn.read().replace( '\n', ' ' )

pyrosetta.init( ''.join( flags ) ) 

ligand_params = pyrosetta.Vector1( [ 'input_files/pNPG.params' ] )
new_res_set = pyrosetta.generate_nonstandard_residue_set( ligand_params )

p = pyrosetta.Pose()
pyrosetta.pose_from_file( p, new_res_set, 'input_files/bglb.pdb' ) 
scorefxn = pyrosetta.create_score_function( 'beta_cst' ) 

add_cst = rosetta.protocols.enzdes.AddOrRemoveMatchCsts()
add_cst.cstfile( 'input_files/pNPG.enzdes.cst' ) 
add_cst.set_cst_action( rosetta.protocols.enzdes.CstAction.ADD_NEW )
add_cst.apply( p ) 

target = int( mutant_name[ 1:-1 ] )
new_res = mutant_name[ -1 ] 
mut = rosetta.protocols.simple_moves.MutateResidue( target, 'LYS' )
mut.apply( p ) 

#pack_task = pyrosetta.standard_packer_task( p )
#pack_task.restrict_to_repacking()

repack = rosetta.protocols.enzdes.EnzRepackMinimize()
repack.set_scorefxn_repack( scorefxn )
repack.set_scorefxn_minimize( scorefxn )
# set all the rest of the fucking options here 

repack.apply( p ) 
p.dump_pdb( 'output_files/{}.pdb'.format( mutant_name) ) 

