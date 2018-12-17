for E in range(300,1600,100):

    lines = '''#!/bin/bash
#PBS -N Lippmann-Schwinger-inf-vol
#PBS -q red
#PBS -l walltime=48:00:00
#PBS -l cput=256:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -o /latticeQCD/raid6/darvish/ls-solver/work/jobs/E_%d.log
#PBS -e /latticeQCD/raid6/darvish/ls-solver/work/jobs/E_%d_err.log


''' % (E,E)
    lines += '''export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export PATH=$PATH:/usr/local/bin

'''
    lines += 'cd /latticeQCD/raid6/darvish/ls-solver/work && /latticeQCD/raid6/darvish/ls-solver/work/test-1b2c %f' % E
    lines += '\n'

    f = open('E_%d_job' % E, 'w')
    f.write(lines)
    f.close()