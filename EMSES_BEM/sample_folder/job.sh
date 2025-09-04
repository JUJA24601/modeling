#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr20001a
#SBATCH -t 72:00:00
#SBATCH --rsc p=114:t=1:c=4:m=4280M
#SBATCH -o %x.%j.out

#============ Shell Script ============
set -x
export LD_LIBRARY_PATH="/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib:/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib:$LD_LIBRARY_PATH"

date
srun -n114 -c4 -l --multi-prog multi.conf
echo ...done
/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python3 generate_xdmf3.py nd1p00_0000.h5 nd2p00_0000.h5 phisp00_0000.h5
/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python3 generate_xdmf3.py ex00_0000.h5 ey00_0000.h5 ez00_0000.h5
/LARGE2/gr10451/share/y-miyake/anaconda3/bin/python3 generate_xdmf3.py j1x00_0000.h5 j1y00_0000.h5 j1z00_0000.h5 j2x00_0000.h5 j2y00_0000.h5 j2z00_0000.h5
date
