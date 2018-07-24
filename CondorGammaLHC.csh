#!/bin/csh 
                                                                                                                                            

#-------------------                                                                                                                                 
# Arguments                                                                                                               
#-------------------                                                                                                                                  

@ p = $1

#-------------------                                                                                                                                
# Variable Defs                                                                                                                                      
#-------------------                                                                                                                                  

set OUT_FILE="/sphenix/user/vassalli/LHCGamma"

set SCRATCH_AREA="$_CONDOR_SCRATCH_DIR"                                                                                                              
#set SCRATCH_AREA="/phenix/scratch/chase"

set SOURCE_PHOTONMAKER="/direct/phenix+u/vassalli/jamieAsym/DirectphotonDijet"

#-------------------                                                                                                                                
# Export Libraries                                                                                                                                   
#-------------------                                                                                                                                  

source /phenix/u/vassalli/.cshrc

#-------------------                                                                                                                                 
# Set Scratch Area                                                                                                                                   
#-------------------                                                                                                                                  

mkdir $SCRATCH_AREA/jamieasym
cp  $SOURCE_PHOTONMAKER $SCRATCH_AREA/jamieasym/

#-------------------                                                                                                                                
# Run Executable  
#-------------------                                                                                                                                  

cd $SCRATCH_AREA/jamieasym
./DirectphotonDijet LHCGamma${1} 95 100 50000
cp LHCGamma${1}* $OUT_FILE


rm -r $SCRATCH_AREA/jamieasym


exit 0
